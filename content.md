# Content

## 1. LinkedIn Post Content

I built an **Interplanetary Mission Studio - Lunar Gateway Planner**, a Streamlit-based mission analysis tool for studying conceptual resupply trajectories from low Earth orbit to the Lunar Gateway region.

The problem statement was simple but technically rich: how can we estimate launch opportunities, transfer trajectories, and delta-v budgets for a Lunar Gateway resupply mission without pretending to build a full operational flight dynamics system?

My solution combines **Poliastro**, **Astropy**, **NumPy**, **Pandas**, **Matplotlib**, and **Streamlit**. The backend solves Earth-to-Moon Lambert transfers, compares a direct transfer against a two-leg Earth -> EM-L1 -> Moon chain, and scans launch windows across departure date and time of flight. I also upgraded the original model by adding **LEO phase optimization**, so the spacecraft no longer departs from one arbitrary point in its parking orbit. The app samples possible parking-orbit phases and selects the lowest injection delta-v.

The data source is primarily **Astropy solar-system ephemerides**, using JPL ephemeris support for Earth-Moon geometry. The app outputs total delta-v, TLI burn, lunar capture burn, Moon arrival v-infinity, best LEO phase angle, injection speed, Earth C3-like energy, ranked launch windows, and a porkchop-style heatmap.

The caveat is important: this is an intermediate-fidelity conceptual planner, not a real mission planner. Gateway is modeled as a circular lunar proxy orbit, not NRHO. The L1 leg is geometric, not CR3BP manifold-based. Solar perturbations, launch-site constraints, finite burns, station-keeping, and vehicle performance are not included.

Still, it is a meaningful step beyond a basic demo: it gives structured mission-design intuition with transparent assumptions.

## 2. Script For My Reel

In this project, I built a **Lunar Gateway Resupply Mission Planner** inside Streamlit.

The problem statement was: given a spacecraft starting from low Earth orbit, can we explore possible transfer options to the Lunar Gateway region, estimate the delta-v budget, and compare launch windows in an interactive way?

I did not want this to be just a static calculation. So I created an app with three main workflows.

First, a **Direct Transfer** tab. This solves an Earth-to-Moon Lambert transfer and estimates the TLI burn, lunar capture burn, total delta-v, Moon arrival v-infinity, injection speed, and C3-like energy.

Second, a **Chain via L1** tab. This models a two-leg path: low Earth orbit to Earth-Moon L1, then L1 to the Moon. It estimates the first departure burn, the matching burn at L1, and the final lunar capture burn.

Third, a **Launch Window Scan** tab. This scans a grid of launch dates and times of flight, ranks the lowest delta-v opportunities, and visualizes the result as a heatmap.

For the data source, I used **Astropy ephemerides with JPL support** to get Earth-Moon geometry, and **Poliastro** for orbital mechanics and Lambert solving.

The key upgrade is **LEO phase optimization**. Instead of assuming one arbitrary point in the parking orbit, the app samples multiple LEO phase angles and chooses the best departure geometry.

The caveat: this is not flight software. Gateway is approximated as a circular lunar orbit, not a true NRHO. L1 is geometric, not a full three-body transfer. But as an intermediate mission-analysis tool, it gives useful and transparent trajectory insight.
