from datetime import date, datetime, time

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from astropy import units as u

from lunar_gateway_mission_planner import (
    compute_lunar_gateway_chain_transfer,
    compute_lunar_gateway_transfer,
    init_ephemeris,
    scan_launch_window,
)


init_ephemeris("jpl")

st.set_page_config(
    page_title="Interplanetary Mission Studio - Lunar Gateway Planner",
    layout="wide",
)


def dt_to_iso(d: date, t: time) -> str:
    dt = datetime(
        year=d.year,
        month=d.month,
        day=d.day,
        hour=t.hour,
        minute=t.minute,
        second=t.second,
    )
    return dt.isoformat()


def qvalue(quantity, unit, digits: int = 3) -> str:
    return f"{quantity.to(unit).value:.{digits}f}"


st.title("Interplanetary Mission Studio")
st.subheader("Lunar Gateway Resupply Mission Planner")

st.sidebar.header("Mission Parameters")

launch_date = st.sidebar.date_input("Launch date", value=date(2030, 1, 10))
launch_time = st.sidebar.time_input("Launch time (UTC)", value=time(12, 0, 0))

alt_leo_km = st.sidebar.number_input(
    "LEO altitude [km]",
    min_value=160.0,
    max_value=2000.0,
    value=400.0,
    step=10.0,
)

gateway_alt_km = st.sidebar.number_input(
    "Gateway proxy orbit altitude [km]",
    min_value=1000.0,
    max_value=10000.0,
    value=3000.0,
    step=100.0,
)

leo_inc_deg = st.sidebar.number_input(
    "LEO inclination [deg]",
    min_value=0.0,
    max_value=98.0,
    value=28.5,
    step=0.5,
)

n_phase_samples = st.sidebar.slider(
    "LEO phase samples",
    min_value=12,
    max_value=144,
    value=36,
    step=12,
    help="Samples positions around the parking orbit and selects the lowest injection delta-v.",
)

st.sidebar.caption(
    "Fidelity: patched-conic Lambert model with JPL ephemerides and parking-orbit "
    "phase optimization. Gateway is still represented by a circular lunar proxy orbit."
)

iso_launch = dt_to_iso(launch_date, launch_time)

tab1, tab2, tab3 = st.tabs(
    ["Direct Transfer", "Chain via L1", "Launch Window Scan"]
)

with tab1:
    st.header("Direct Earth -> Moon Transfer")

    tof_days = st.number_input(
        "Time of flight [days]",
        min_value=2.0,
        max_value=10.0,
        value=3.5,
        step=0.25,
    )

    if st.button("Compute direct transfer", key="direct_btn"):
        with st.spinner("Solving Lambert transfers across the LEO phase grid..."):
            try:
                res = compute_lunar_gateway_transfer(
                    launch_date=iso_launch,
                    tof_days=tof_days,
                    alt_leo_km=alt_leo_km,
                    gateway_alt_km=gateway_alt_km,
                    leo_inc_deg=leo_inc_deg,
                    n_phase_samples=n_phase_samples,
                )

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Total delta-v", f"{qvalue(res['dv_total'], u.km / u.s)} km/s")
                m2.metric("TLI burn", f"{qvalue(res['dv_TLI'], u.km / u.s)} km/s")
                m3.metric("Lunar capture", f"{qvalue(res['dv_LOI_gateway'], u.km / u.s)} km/s")
                m4.metric("Moon v_inf", f"{qvalue(res['v_inf_moon'], u.km / u.s)} km/s")

                col1, col2, col3 = st.columns(3)
                with col1:
                    st.markdown("### Epochs")
                    st.write(f"Launch epoch: `{res['launch_epoch'].iso}`")
                    st.write(f"Arrival epoch: `{res['arrival_epoch'].iso}`")
                    st.write(f"Time of flight: `{res['tof'].to(u.day):.3f}`")

                with col2:
                    st.markdown("### Departure Geometry")
                    st.write(f"Best LEO phase: **{res['leo_phase_deg']:.1f} deg**")
                    st.write(f"Phase samples searched: **{res['n_phase_samples']}**")
                    st.write(
                        "Parking speed: "
                        f"**{qvalue(res['parking_speed'], u.km / u.s)} km/s**"
                    )
                    st.write(
                        "Injection speed: "
                        f"**{qvalue(res['injection_speed'], u.km / u.s)} km/s**"
                    )

                with col3:
                    st.markdown("### Energy")
                    st.write(
                        "Earth C3-like energy: "
                        f"**{qvalue(res['earth_c3'], u.km**2 / u.s**2)} km2/s2**"
                    )
                    st.write(
                        "This is a geocentric energy diagnostic. Negative values "
                        "can occur for bound Earth-centered lunar transfer arcs."
                    )

            except Exception as e:
                st.error(f"Failed to solve transfer: {e}")

with tab2:
    st.header("Two-Leg Chain: Earth -> EM-L1 -> Lunar Gateway")

    col_a, col_b = st.columns(2)
    with col_a:
        tof1 = st.number_input(
            "Leg 1: LEO -> L1 TOF [days]",
            min_value=0.5,
            max_value=5.0,
            value=1.5,
            step=0.25,
        )
    with col_b:
        tof2 = st.number_input(
            "Leg 2: L1 -> Moon TOF [days]",
            min_value=1.0,
            max_value=7.0,
            value=2.0,
            step=0.25,
        )

    if st.button("Compute chain transfer", key="chain_btn"):
        with st.spinner("Solving two Lambert legs with optimized LEO departure phase..."):
            try:
                res = compute_lunar_gateway_chain_transfer(
                    launch_date=iso_launch,
                    tof_leg1_days=tof1,
                    tof_leg2_days=tof2,
                    alt_leo_km=alt_leo_km,
                    gateway_alt_km=gateway_alt_km,
                    leo_inc_deg=leo_inc_deg,
                    n_phase_samples=n_phase_samples,
                )

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Total delta-v", f"{qvalue(res['dv_total'], u.km / u.s)} km/s")
                m2.metric("LEO -> L1 burn", f"{qvalue(res['leg1']['dv_TLI'], u.km / u.s)} km/s")
                m3.metric("L1 match burn", f"{qvalue(res['leg2']['dv_mid'], u.km / u.s)} km/s")
                m4.metric("Lunar capture", f"{qvalue(res['dv_LOI_gateway'], u.km / u.s)} km/s")

                col1, col2, col3 = st.columns(3)
                with col1:
                    st.markdown("### Epochs")
                    st.write(f"Launch epoch: `{res['launch_epoch'].iso}`")
                    st.write(f"L1 epoch: `{res['L1_epoch'].iso}`")
                    st.write(f"Arrival epoch: `{res['arrival_epoch'].iso}`")
                    st.write(f"Leg 1 TOF: `{res['leg1']['tof'].to(u.day):.3f}`")
                    st.write(f"Leg 2 TOF: `{res['leg2']['tof'].to(u.day):.3f}`")

                with col2:
                    st.markdown("### Departure Geometry")
                    st.write(f"Best LEO phase: **{res['leg1']['leo_phase_deg']:.1f} deg**")
                    st.write(f"Phase samples searched: **{res['n_phase_samples']}**")
                    st.write(
                        "Injection speed: "
                        f"**{qvalue(res['leg1']['injection_speed'], u.km / u.s)} km/s**"
                    )
                    st.write(
                        "Earth C3-like energy: "
                        f"**{qvalue(res['leg1']['earth_c3'], u.km**2 / u.s**2)} km2/s2**"
                    )

                with col3:
                    st.markdown("### Arrival")
                    st.write(
                        "Moon v_inf on leg 2: "
                        f"**{qvalue(res['leg2']['v_inf_moon'], u.km / u.s)} km/s**"
                    )
                    st.write(
                        "L1 remains a geometric waypoint in this model, not a full "
                        "three-body invariant-manifold transfer."
                    )

            except Exception as e:
                st.error(f"Failed to solve chain transfer: {e}")

with tab3:
    st.header("Launch Window Scan")

    st.caption(
        "Grid search over launch dates and flight times. Each grid cell includes "
        "a parking-orbit phase search, so large grids can take time."
    )

    col1, col2 = st.columns(2)
    with col1:
        lw_start = st.date_input("Scan start date", value=date(2030, 1, 1))
        lw_end = st.date_input("Scan end date", value=date(2030, 1, 15))

    with col2:
        tof_min = st.number_input("Minimum TOF [days]", 2.0, 10.0, 3.0, 0.25)
        tof_max = st.number_input("Maximum TOF [days]", 2.0, 15.0, 5.0, 0.25)
        n_launch = st.slider("Number of launch samples", 5, 40, 12, 1)
        n_tof = st.slider("Number of TOF samples", 5, 40, 10, 1)

    if st.button("Run launch window scan", key="scan_btn"):
        if lw_end <= lw_start:
            st.error("Scan end date must be after start date.")
        elif tof_max <= tof_min:
            st.error("Maximum TOF must be greater than minimum TOF.")
        else:
            with st.spinner("Scanning launch window and optimizing LEO phase..."):
                try:
                    df = scan_launch_window(
                        launch_start=lw_start.isoformat(),
                        launch_end=lw_end.isoformat(),
                        tof_min_days=tof_min,
                        tof_max_days=tof_max,
                        n_launch=n_launch,
                        n_tof=n_tof,
                        alt_leo_km=alt_leo_km,
                        gateway_alt_km=gateway_alt_km,
                        leo_inc_deg=leo_inc_deg,
                        n_phase_samples=n_phase_samples,
                    )

                    if df.empty:
                        st.warning("No valid Lambert solutions found for this window.")
                    else:
                        st.subheader("Top 10 lowest-delta-v transfers")
                        st.dataframe(df.head(10), use_container_width=True)

                        df_heat = df.copy()
                        df_heat["launch_day"] = pd.to_datetime(
                            df_heat["launch_epoch"]
                        ).dt.date
                        df_heat = (
                            df_heat.sort_values("dv_total_kms")
                            .groupby(["launch_day", "tof_days"], as_index=False)
                            .first()
                        )
                        pivot = df_heat.pivot(
                            index="tof_days",
                            columns="launch_day",
                            values="dv_total_kms",
                        )

                        fig, ax = plt.subplots(figsize=(10, 5))
                        im = ax.imshow(pivot.values, aspect="auto", origin="lower")
                        ax.set_xticks(range(pivot.shape[1]))
                        ax.set_xticklabels(
                            pivot.columns.astype(str), rotation=45, ha="right"
                        )
                        ax.set_yticks(range(pivot.shape[0]))
                        ax.set_yticklabels([f"{v:.2f}" for v in pivot.index])
                        ax.set_xlabel("Launch date")
                        ax.set_ylabel("TOF [days]")
                        ax.set_title("Total delta-v [km/s] - direct transfers")

                        cbar = fig.colorbar(im, ax=ax)
                        cbar.set_label("Total delta-v [km/s]")
                        st.pyplot(fig)

                except Exception as e:
                    st.error(f"Scan failed: {e}")
