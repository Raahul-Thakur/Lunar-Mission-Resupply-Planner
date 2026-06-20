"""
Mission Planner for a Lunar Gateway Resupply Mission
----------------------------------------------------

Now supports:
- Direct Earth → Moon transfers
- Two-leg chain: Earth → EM-L1 → Moon (Gateway capture)

Dependencies:
    pip install poliastro astropy numpy pandas matplotlib
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris

from poliastro.bodies import Earth, Moon
from poliastro.ephem import Ephem
from poliastro.twobody import Orbit
from poliastro.iod.izzo import lambert
from poliastro.plotting.porkchop import PorkchopPlotter
from poliastro.util import time_range
import matplotlib.pyplot as plt


def _norm(vector: u.Quantity) -> u.Quantity:
    """Return Euclidean norm while preserving Astropy units."""
    return np.linalg.norm(vector.to(vector.unit).value) * vector.unit


def _earth_c3(v_geo: u.Quantity, r_geo: u.Quantity) -> u.Quantity:
    """
    Twice the geocentric specific orbital energy for the departure arc.

    Positive values are classical C3. Lunar transfers can also produce
    negative values because the Earth-centered conic is not necessarily
    hyperbolic.
    """
    speed = v_geo.to(u.km / u.s)
    speed2 = np.dot(speed.value, speed.value) * (u.km**2 / u.s**2)
    radius = _norm(r_geo.to(u.km))
    return (speed2 - 2 * Earth.k.to(u.km**3 / u.s**2) / radius).to(
        u.km**2 / u.s**2
    )


# ----------------------------------------------------------------------
# 0. Ephemeris init
# ----------------------------------------------------------------------

def init_ephemeris(ephem: str = "jpl") -> None:
    """
    Set the solar system ephemeris backend.
    'jpl' is more accurate but requires an extra download the first time.
    """
    solar_system_ephemeris.set(ephem)


# ----------------------------------------------------------------------
# 1. Basic orbits and L1 approximation
# ----------------------------------------------------------------------

def leo_parking_orbit(
    epoch: Time,
    alt_leo: float = 400.0,
    inc: float = 28.5,
    arglat: u.Quantity = 0.0 * u.deg,
) -> Orbit:
    """
    Create a circular LEO parking orbit.

    Parameters
    ----------
    epoch : astropy.time.Time
        Departure epoch.
    alt_leo : float
        Altitude above Earth's surface [km].
    inc : float
        Inclination [deg].

    Returns
    -------
    Orbit
        LEO orbit object around Earth.
    """
    return Orbit.circular(
        Earth,
        alt_leo * u.km,
        inc=inc * u.deg,
        arglat=arglat,
        epoch=epoch,
    )


def moon_geocentric_orbit(epoch: Time) -> Orbit:
    """
    Get the Moon's osculating orbit around Earth at 'epoch'.
    """
    moon_ephem = Ephem.from_body(Moon, epoch, attractor=Earth)
    return Orbit.from_ephem(Earth, moon_ephem, epoch)


def earth_moon_L1_position(epoch: Time) -> u.Quantity:
    """
    Approximate the Earth–Moon L1 point position in Earth-centered frame.

    Theory:
        L1 lies on the Earth–Moon line, at a distance ~ d_EM - d_L1 from Earth,
        where d_L1 ≈ a * (μ/3)^(1/3), μ = m_moon/(m_earth + m_moon).

    We use actual Earth→Moon vector length at 'epoch' and the mass ratio from
    poliastro bodies.

    Returns
    -------
    r_L1 : (3,) Quantity
        Position vector of L1 in km in the same frame as Moon's geocentric orbit.
    """
    moon_orbit = moon_geocentric_orbit(epoch)
    r_moon = moon_orbit.r.to(u.km)  # Earth→Moon vector

    # Mass ratio μ = m_moon / (m_earth + m_moon)
    mu_mass = Moon.mass / (Earth.mass + Moon.mass)

    # Earth–Moon distance (scalar)
    d_em = np.linalg.norm(r_moon.value) * u.km

    # Distance from Moon toward Earth to L1
    # classical collinear L1 approximation: d_L1 ≈ a * (μ/3)^(1/3)
    d_L1 = d_em * (mu_mass / 3.0) ** (1.0 / 3.0)

    # Unit vector from Earth to Moon
    r_hat_em = r_moon / d_em

    # L1 lies between Earth and Moon:
    # position measured from Earth: r_L1 = r_moon - d_L1 * r_hat_em
    r_L1 = r_moon - d_L1 * r_hat_em
    return r_L1.to(u.km)


def best_leo_departure_lambert(
    target_r: u.Quantity,
    tof: u.Quantity,
    epoch: Time,
    alt_leo_km: float = 400.0,
    leo_inc_deg: float = 28.5,
    n_phase_samples: int = 72,
) -> dict:
    """
    Search the circular parking orbit for the lowest-delta-v Lambert departure.

    This keeps the model lightweight but avoids using one arbitrary LEO point
    for every launch opportunity.
    """
    if n_phase_samples < 1:
        raise ValueError("n_phase_samples must be at least 1")

    best = None
    phases = np.linspace(-180.0, 180.0, n_phase_samples, endpoint=False)

    for phase_deg in phases:
        leo = leo_parking_orbit(
            epoch,
            alt_leo=alt_leo_km,
            inc=leo_inc_deg,
            arglat=phase_deg * u.deg,
        )
        v_depart, v_arrive = lambert(Earth.k, leo.r, target_r, tof)
        dv_depart = _norm((v_depart - leo.v).to(u.km / u.s))

        if best is None or dv_depart < best["dv_depart"]:
            best = {
                "leo": leo,
                "phase_deg": phase_deg % 360.0,
                "v_depart": v_depart,
                "v_arrive": v_arrive,
                "dv_depart": dv_depart,
                "parking_speed": _norm(leo.v.to(u.km / u.s)),
                "injection_speed": _norm(v_depart.to(u.km / u.s)),
                "earth_c3": _earth_c3(v_depart, leo.r),
            }

    return best


# ----------------------------------------------------------------------
# 2. Direct Earth → Moon transfer + Δv budget
# ----------------------------------------------------------------------

def compute_lunar_gateway_transfer(
    launch_date: str,
    tof_days: float,
    alt_leo_km: float = 400.0,
    gateway_alt_km: float = 3000.0,
    leo_inc_deg: float = 28.5,
    n_phase_samples: int = 72,
) -> dict:
    """
    Compute a single Earth→Moon transfer and approximate Δv budget
    for a Lunar Gateway resupply mission (direct transfer).

    Returns
    -------
    dict with keys:
        launch_epoch, arrival_epoch, tof,
        dv_TLI, dv_LOI_gateway, dv_total, v_inf_moon
    """
    # Epochs
    launch_epoch = Time(launch_date, scale="tdb")
    tof = tof_days * u.day
    arrival_epoch = launch_epoch + tof

    # Initial and target
    moon_orbit_arr = moon_geocentric_orbit(arrival_epoch)
    r2 = moon_orbit_arr.r

    departure = best_leo_departure_lambert(
        target_r=r2,
        tof=tof,
        epoch=launch_epoch,
        alt_leo_km=alt_leo_km,
        leo_inc_deg=leo_inc_deg,
        n_phase_samples=n_phase_samples,
    )

    # Lambert around Earth
    v2_trans = departure["v_arrive"]

    # Injection Δv
    dv_TLI = departure["dv_depart"]

    # Hyperbolic excess at Moon
    v_inf_vec = v2_trans - moon_orbit_arr.v
    v_inf = _norm(v_inf_vec.to(u.km / u.s))

    # Gateway orbit around Moon
    mu_moon = Moon.k.to(u.km**3 / u.s**2)
    r_gateway = (Moon.R + gateway_alt_km * u.km).to(u.km)

    v_peri_hyp = np.sqrt(v_inf**2 + 2 * mu_moon / r_gateway)
    v_circ_gateway = np.sqrt(mu_moon / r_gateway)

    dv_LOI_gateway = (v_peri_hyp - v_circ_gateway).to(u.km / u.s)
    dv_total = dv_TLI + dv_LOI_gateway

    return {
        "launch_epoch": launch_epoch,
        "arrival_epoch": arrival_epoch,
        "tof": tof,
        "dv_TLI": dv_TLI,
        "dv_LOI_gateway": dv_LOI_gateway,
        "dv_total": dv_total,
        "v_inf_moon": v_inf,
        "leo_phase_deg": departure["phase_deg"],
        "parking_speed": departure["parking_speed"],
        "injection_speed": departure["injection_speed"],
        "earth_c3": departure["earth_c3"],
        "n_phase_samples": n_phase_samples,
    }


# ----------------------------------------------------------------------
# 3. Two-leg chain: Earth → L1 → Moon (Gateway capture)
# ----------------------------------------------------------------------

def compute_lunar_gateway_chain_transfer(
    launch_date: str,
    tof_leg1_days: float,
    tof_leg2_days: float,
    alt_leo_km: float = 400.0,
    gateway_alt_km: float = 3000.0,
    leo_inc_deg: float = 28.5,
    n_phase_samples: int = 72,
) -> dict:
    """
    Two-leg chain transfer:
        Leg 1: LEO (Earth) → EM-L1
        Leg 2: EM-L1 → Moon
        Final: capture into Lunar Gateway orbit.

    Each leg is solved via Lambert in the Earth-centered frame.
    The mid-course maneuver Δv at L1 matches the two transfer arcs.

    Parameters
    ----------
    launch_date : str
        Launch date (ISO).
    tof_leg1_days : float
        Time of flight for LEO → L1 [days].
    tof_leg2_days : float
        Time of flight for L1 → Moon [days].

    Returns
    -------
    dict
        {
          "leg1": {...},
          "leg2": {...},
          "dv_LOI_gateway": Quantity,
          "dv_total": Quantity,
          "launch_epoch": Time,
          "L1_epoch": Time,
          "arrival_epoch": Time,
        }
    """
    launch_epoch = Time(launch_date, scale="tdb")
    tof1 = tof_leg1_days * u.day
    L1_epoch = launch_epoch + tof1
    tof2 = tof_leg2_days * u.day
    arrival_epoch = L1_epoch + tof2

    # States
    r_L1 = earth_moon_L1_position(L1_epoch)
    moon_orbit_arr = moon_geocentric_orbit(arrival_epoch)
    r_moon_arr = moon_orbit_arr.r

    # Leg 1: Lambert LEO → L1 (Earth centered)
    leg1_departure = best_leo_departure_lambert(
        target_r=r_L1,
        tof=tof1,
        epoch=launch_epoch,
        alt_leo_km=alt_leo_km,
        leo_inc_deg=leo_inc_deg,
        n_phase_samples=n_phase_samples,
    )
    v1_leg1 = leg1_departure["v_depart"]
    v2_leg1 = leg1_departure["v_arrive"]

    dv1_TLI = leg1_departure["dv_depart"]

    # Leg 2: Lambert L1 → Moon (Earth centered)
    v1_leg2, v2_leg2 = lambert(Earth.k, r_L1, r_moon_arr, tof2)

    # Mid-course Δv at L1 (match end of leg1 to start of leg2)
    dv_mid_vec = v1_leg2 - v2_leg1
    dv_mid = _norm(dv_mid_vec.to(u.km / u.s))

    # Relative velocity at Moon arrival for capture
    v_inf_vec = v2_leg2 - moon_orbit_arr.v
    v_inf = _norm(v_inf_vec.to(u.km / u.s))

    # Lunar Gateway capture (same as direct case)
    mu_moon = Moon.k.to(u.km**3 / u.s**2)
    r_gateway = (Moon.R + gateway_alt_km * u.km).to(u.km)

    v_peri_hyp = np.sqrt(v_inf**2 + 2 * mu_moon / r_gateway)
    v_circ_gateway = np.sqrt(mu_moon / r_gateway)

    dv_LOI_gateway = (v_peri_hyp - v_circ_gateway).to(u.km / u.s)

    dv_total = dv1_TLI + dv_mid + dv_LOI_gateway

    return {
        "launch_epoch": launch_epoch,
        "L1_epoch": L1_epoch,
        "arrival_epoch": arrival_epoch,
        "leg1": {
            "tof": tof1,
            "dv_TLI": dv1_TLI,
            "v_depart": v1_leg1,
            "v_arrive_L1": v2_leg1,
            "leo_phase_deg": leg1_departure["phase_deg"],
            "parking_speed": leg1_departure["parking_speed"],
            "injection_speed": leg1_departure["injection_speed"],
            "earth_c3": leg1_departure["earth_c3"],
        },
        "leg2": {
            "tof": tof2,
            "dv_mid": dv_mid,
            "v_depart_L1": v1_leg2,
            "v_arrive_moon": v2_leg2,
            "v_inf_moon": v_inf,
        },
        "dv_LOI_gateway": dv_LOI_gateway,
        "dv_total": dv_total,
        "n_phase_samples": n_phase_samples,
    }


# ----------------------------------------------------------------------
# 4. Launch window scan (direct transfers, Δv grid)
# ----------------------------------------------------------------------

def scan_launch_window(
    launch_start: str,
    launch_end: str,
    tof_min_days: float = 3.0,
    tof_max_days: float = 7.0,
    n_launch: int = 15,
    n_tof: int = 15,
    alt_leo_km: float = 400.0,
    gateway_alt_km: float = 3000.0,
    leo_inc_deg: float = 28.5,
    n_phase_samples: int = 36,
) -> pd.DataFrame:
    """
    Scan a grid of [launch date, time-of-flight] and compute Δv budgets
    for direct Earth→Moon transfers.

    Acts like a simple porkchop in Δv space.
    """
    launch_span = Time(launch_start) + np.linspace(
        0.0,
        (Time(launch_end) - Time(launch_start)).to(u.day).value,
        n_launch,
    ) * u.day

    tof_span = np.linspace(tof_min_days, tof_max_days, n_tof)

    rows = []

    for t_launch in launch_span:
        for tof_days in tof_span:
            try:
                result = compute_lunar_gateway_transfer(
                    launch_date=t_launch.iso,
                    tof_days=tof_days,
                    alt_leo_km=alt_leo_km,
                    gateway_alt_km=gateway_alt_km,
                    leo_inc_deg=leo_inc_deg,
                    n_phase_samples=n_phase_samples,
                )
                rows.append(
                    {
                        "launch_epoch": result["launch_epoch"].iso,
                        "arrival_epoch": result["arrival_epoch"].iso,
                        "tof_days": result["tof"].to(u.day).value,
                        "dv_TLI_kms": result["dv_TLI"].to(u.km / u.s).value,
                        "dv_LOI_gateway_kms": result[
                            "dv_LOI_gateway"
                        ].to(u.km / u.s).value,
                        "dv_total_kms": result["dv_total"].to(u.km / u.s).value,
                        "v_inf_moon_kms": result["v_inf_moon"].to(u.km / u.s).value,
                        "earth_c3_km2_s2": result["earth_c3"].to(
                            u.km**2 / u.s**2
                        ).value,
                        "leo_phase_deg": result["leo_phase_deg"],
                    }
                )
            except Exception as e:
                print(f"Warning: Lambert failed at {t_launch.iso}, TOF={tof_days} d: {e}")

    df = pd.DataFrame(rows)
    if not df.empty:
        df.sort_values("dv_total_kms", inplace=True, ignore_index=True)
    return df


# ----------------------------------------------------------------------
# 5. Optional: classical v_inf / C3 porkchop
# ----------------------------------------------------------------------

def plot_earth_moon_porkchop(
    launch_start: str,
    launch_end: str,
    arrival_start: str,
    arrival_end: str,
    n_launch: int = 40,
    n_arrival: int = 40,
) -> None:
    """
    Use poliastro's PorkchopPlotter to visualize Earth→Moon launch windows.
    """
    launch_span = time_range(
        launch_start,
        end=launch_end,
        periods=n_launch,
    )
    arrival_span = time_range(
        arrival_start,
        end=arrival_end,
        periods=n_arrival,
    )

    pork = PorkchopPlotter(Earth, Moon, launch_span, arrival_span)
    fig, ax = pork.plot()
    ax.set_title("Earth → Moon Porkchop Plot (C3 / v_inf contours)")
    plt.show()


# ----------------------------------------------------------------------
# 6. CLI test
# ----------------------------------------------------------------------

if __name__ == "__main__":
    init_ephemeris("jpl")

    print("=== Direct Transfer Example ===")
    direct = compute_lunar_gateway_transfer(
        launch_date="2030-01-10 12:00:00",
        tof_days=3.5,
    )
    print(f"Δv_total (direct) : {direct['dv_total']:.3f}")

    print("\n=== Chain Transfer Example (Earth → L1 → Moon) ===")
    chain = compute_lunar_gateway_chain_transfer(
        launch_date="2030-01-10 12:00:00",
        tof_leg1_days=1.5,
        tof_leg2_days=2.0,
    )
    print(f"Δv_total (chain)  : {chain['dv_total']:.3f}")
    print(f"  - dv_TLI (leg1) : {chain['leg1']['dv_TLI']:.3f}")
    print(f"  - dv_mid (L1)   : {chain['leg2']['dv_mid']:.3f}")
    print(f"  - dv_LOI        : {chain['dv_LOI_gateway']:.3f}")
