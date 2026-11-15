"""
Converts Perseverance (Mars 2020) traverse + waypoints (NASA/JPL MMGIS GeoJSON)
into small files for the web app, and places sample sites at the rover's
position on the sol each sample tube was sealed.
Usage: python rover.py <M20_waypoints.json> <M20_traverse.json> <out_dir>
"""
import sys, json
wp_path, tr_path, out = sys.argv[1:4]
wp = [f['properties'] for f in json.load(open(wp_path))['features']]
tr = json.load(open(tr_path))['features']

# 1. Traverse: one simplified line (lon, lat) for the whole mission
line = []
for f in tr:
    for x, y, *_ in f['geometry']['coordinates']:
        p = [round(x, 5), round(y, 5)]
        # keep a point only if it is >~10 m (0.0002 deg) from the previous kept one
        if not line or abs(p[0]-line[-1][0]) + abs(p[1]-line[-1][1]) > 0.0002:
            line.append(p)
# 2. Waypoints: sol + position + distance driven so far
pts = [[round(p['lon'], 5), round(p['lat'], 5), p['sol'], p['dist_km']] for p in wp]

# 3. Sample sites (NASA "Mars Rock Samples" page: number, name, sol sealed, rock type)
SAMPLES = [
  (9,  'Atsá',            377,  'Sid',              'Igneous',   'Crater floor lava rock: dates the floor and records its interaction with water.'),
  (11, 'Skyland',         495,  'Skinner Ridge',    'Sedimentary','Delta front sandstone. Its twin now sits in the Three Forks sample depot.'),
  (12, 'Hazeltop',        509,  'Wildcat Ridge',    'Sedimentary','Fine-grained mudstone from an ancient lake: best kind of rock to preserve organics.'),
  (13, 'Bearwallow',      575,  'Wildcat Ridge',    'Sedimentary','Paired lake-mudstone core from the delta front.'),
  (15, 'Mageik',          619,  'Amalik',           'Sedimentary','Sulfate-bearing delta rock: evidence of evaporating water.'),
  (17, 'Atmo Mountain',   634,  'Observation Mountain','Regolith','Loose surface material: key for dust hazards to astronauts and equipment.'),
  (20, 'Otis Peak',       832,  'Emerald Lake',     'Sedimentary','Phosphate-bearing rock from the upper delta fan.'),
  (21, 'Pilot Mountain',  913,  'Dream Lake',       'Sedimentary','Coarse river deposit carrying pebbles from outside Jezero.'),
  (22, 'Pelican Point',   923,  'Hans Amundsen Memorial Workspace','Sedimentary','Upper fan boulder deposit from a high-energy flood.'),
  (24, 'Comet Geyser',    1088, 'Bunsen Peak',      'Silica-cemented carbonate','Carbonate from the crater margin: records past climate and water chemistry.'),
  (25, 'Sapphire Canyon', 1215, 'Cheyava Falls',    'Sedimentary','Neretva Vallis riverbed rock with "leopard spots": potential biosignature (Nature, 2025).'),
  (26, 'Silver Mountain', 1401, 'Shallow Bay',      'Igneous (ancient)','First crater-rim sample: likely the oldest rock collected.'),
]
samples = []
for num, name, sol, feat, rock, why in SAMPLES:
    at = [p for p in wp if p['sol'] <= sol][-1]
    samples.append(dict(num=num, name=name, sol=sol, feature=feat, rock=rock, why=why,
                        lon=round(at['lon'], 5), lat=round(at['lat'], 5), at_sol=at['sol']))
json.dump(dict(traverse=line, waypoints=pts, samples=samples,
               last_sol=wp[-1]['sol'], total_km=wp[-1]['dist_km']),
          open(f'{out}/rover.json', 'w'), separators=(',', ':'))
print(len(line), 'line pts,', len(pts), 'waypoints')
for s in samples: print(s['num'], s['name'], s['sol'], s['at_sol'], s['lon'], s['lat'])
