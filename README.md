# Marswalk Planner · Jezero Crater

**NASA Space Apps Challenge 2026 · "Interplanetary Survival Guide: Martian Map"**
Team *Atoll to Ares* (Zeyk Ahmed, Maldives) · NASA Space Apps Batticaloa

A mission-planning map for a human Marswalk (EVA) at Jezero Crater. Place a habitat, add science
targets, and the planner computes and compares walking routes over real Mars terrain, with safety limits.

**Live demo:** https://<your-username>.github.io/marswalk-planner/

## Features
- Layered map: MRO CTX imagery, CTX elevation, hillshade, slope, hazard classes
- Least-cost routing (A*) habitat → targets → habitat, over a 20 m terrain cost grid
- Route comparison: **Fastest** vs **Safest** vs **Most science** (distance, walk time, climb, max slope, hazardous ground, sample sites passed)
- Perseverance traverse (sol 0–2022, 45 km) and 12 sample sites with science cards; the Most science route auto-adds detours to sample sites while staying within the EVA budget
- Safety: EVA time budget, max distance from habitat, safe-return zone, turn-back point
- Elevation profile coloured by hazard

## Method
1. `scripts/preprocess.py` resamples the CTX DTM and orthoimage onto one 20 m grid
   (equirectangular projection, Mars sphere R = 3,396,190 m, true scale at 18.4663°N).
2. Slope (degrees) and roughness (local relief in a 60 m window) are computed from elevation.
3. Hazard classes: easy < 8°, moderate 8–15°, difficult 15–25° (or top 2% roughness), no-go > 25°.
4. Walking time per step uses Tobler's hiking function × 0.6 (suited-astronaut assumption).
   Route profiles multiply time by a hazard penalty.
5. `scripts/rover.py` converts the Perseverance traverse/waypoints and places each sample at the rover's
   position on the sol its tube was sealed (approximate location).
6. Most science route: greedy insertion of the sample site with the smallest detour, accepted only if
   the walk stays within 90% of the EVA budget and the distance limit.
7. A reverse Dijkstra search from the habitat gives return time from every cell → safe-return
   zone and turn-back point.

Assumptions/limits: thresholds are design assumptions, not NASA EVA requirements; 20 m data
cannot see boulders smaller than ~20 m. Educational prototype, not for real mission use.

## Run locally
```
python scripts/rover.py M20_waypoints.json M20_traverse.json docs/data
python scripts/preprocess.py M20_JezeroCrater_CTXDEM_20m.tif JEZ_ctx_B_soc_008_orthoMosaic_6m_Eqc_latTs0_lon0.tif docs/data
cd docs && python -m http.server
```

## References
### Data
- Calef, F. (2021). *Mars 2020 Science Investigation CTX DEM Mosaic* (20 m), USGS Astrogeology / PDS.
  https://astrogeology.usgs.gov/search/map/mars_2020_science_investigation_ctx_dem_mosaic
- Fergason, R. L., et al. (2020). *Mars 2020 Terrain Relative Navigation CTX Orthorectified Image Mosaic* (6 m). USGS Astrogeology.
  https://astrogeology.usgs.gov/search/map/mars_2020_terrain_relative_navigation_context_camera_orthorectified_image_mosaic
- Fergason, R. L., et al. (2020). *Mars 2020 Terrain Relative Navigation Flight Product Generation*. LPSC 51, abstract 2020.
- Malin, M. C., et al. (2007). Context Camera Investigation on board the Mars Reconnaissance Orbiter. *JGR Planets*.
- NASA/JPL-Caltech, Mars 2020 Perseverance traverse and waypoints (MMGIS "Where is Perseverance?" map), M20_waypoints.json & M20_traverse.json.
- NASA, *Mars Rock Samples* (sample names, numbers, sols, rock types). https://science.nasa.gov/mission/mars-2020-perseverance/mars-rock-samples
- Missions: NASA Mars Reconnaissance Orbiter (CTX), NASA Mars 2020 Perseverance.
### Methods
- Tobler, W. (1993). Three presentations on geographical analysis and modeling (hiking function).
- Hart, Nilsson & Raphael (1968). A* search. Dijkstra (1959).
### Software
- Leaflet 1.9.4 (BSD-2), rasterio, NumPy, SciPy, Pillow, GDAL/PROJ.
- Fonts: Space Grotesk, JetBrains Mono (SIL OFL, Google Fonts).

## AI disclosure
Code was written with the help of Claude (Anthropic), an AI assistant, during the hackathon.
Concept, design decisions, testing, and data selection by Zeyk Ahmed.
