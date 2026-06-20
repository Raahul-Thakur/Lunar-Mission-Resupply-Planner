# Interplanetary Mission Studio - Lunar Gateway Planner

A Streamlit-based mission analysis tool for conceptual Lunar Gateway resupply trajectory studies.

The app uses Poliastro and Astropy to evaluate Earth-to-Moon transfers, estimate delta-v budgets, compare direct and L1-waypoint trajectories, and scan launch windows. It is designed as an intermediate-fidelity engineering study tool: more capable than a classroom demo, but not a flight-dynamics replacement for high-fidelity NRHO or CR3BP mission design.

## Features

- Direct Earth-to-Moon Lambert transfer analysis
- Two-leg Earth -> EM-L1 -> Moon transfer chain
- Launch-window grid scan over departure date and time of flight
- JPL ephemeris support through Astropy
- Circular LEO parking-orbit model with optimized departure phase sampling
- Delta-v breakdown:
  - TLI / departure burn
  - L1 matching burn for the chained trajectory
  - lunar capture burn into a Gateway proxy orbit
- Additional diagnostics:
  - best LEO phase angle
  - parking orbit speed
  - injection speed
  - Earth C3-like energy
  - lunar arrival v-infinity

## Scientific Scope

This project is suitable for conceptual trade studies, early trajectory intuition, and comparing launch windows or transfer options.

It is not a real operational mission planner. The model still simplifies several important effects:

- Gateway is represented as a circular lunar proxy orbit, not a true NRHO.
- Transfers use patched-conic Earth-centered Lambert arcs.
- The L1 waypoint is geometric, not a full CR3BP invariant-manifold trajectory.
- Solar perturbations, station-keeping, launch-site constraints, finite burns, navigation errors, and vehicle performance limits are not modeled.

## Project Structure

```text
app.py
lunar_gateway_mission_planner.py
requirements.txt
README.md
```

- `app.py` contains the Streamlit interface.
- `lunar_gateway_mission_planner.py` contains the mission analysis functions.
- `requirements.txt` pins the environment used for this version.

## Installation

Create and activate a Python 3.10 virtual environment, then install dependencies:

```bash
pip install -r requirements.txt
```

## Run

```bash
streamlit run app.py
```

Then open:

```text
http://localhost:8501
```

The first run may download JPL ephemeris data through Astropy.

## Notes

Astropy may print `dubious year` warnings for future dates such as 2030. These are time-scale precision warnings, not fatal application errors.

Increasing the LEO phase sample count improves the departure geometry search, but also makes each transfer and launch-window scan slower.
