/* Marswalk Planner – routing + UI
 * Grid: 20 m cells, row 0 = north edge. Leaflet runs in EPSG:4326 (plain lon/lat degrees),
 * which matches the equirectangular Mars data, so pixels map linearly to lat/lon. */

const SUIT_FACTOR = 0.6;           // suited astronaut speed vs. Earth hiker (assumption)
const NOGO = 3;                    // hazard class that is never entered

// Route profiles: hazard penalty multipliers for classes [easy, moderate, difficult, no-go]
const PROFILES = [
  { id: 'fast', name: 'Fastest', color: '#ff7a3d', mult: [1, 1.15, 1.6, Infinity] },
  { id: 'safe', name: 'Safest',  color: '#3ed6e0', mult: [1, 3, 12, Infinity] },
  // Most science: safe-style costs, plus automatic detours to nearby Perseverance sample sites
  { id: 'sci',  name: 'Most science', color: '#b98cff', mult: [1, 2, 8, Infinity], science: true },
];
let ROVER, SITES = [];   // Perseverance data; SITES = sample locations grouped by position

let G;                 // grid: {w,h,cell,elev,slope,haz,bounds}
let map, habitat, targets = [], routes = [], selected = 'safe', mode = 'target';
let returnTime = null, zoneLayer = null, routeLayers = [], turnMarker = null;
const state = { budgetH: 4, radiusKm: 5 };

/* ---------- loading ---------- */
async function loadGrid() {
  const meta = await (await fetch('data/meta.json')).json();
  const bin = async (f, T) => new T(await (await fetch('data/' + f)).arrayBuffer());
  G = { w: meta.width, h: meta.height, cell: meta.cell_m, b: meta.bounds,
        elev: await bin('elev.f32', Float32Array), slope: await bin('slope.u16', Uint16Array),
        haz: await bin('haz.u8', Uint8Array) };
  ROVER = await (await fetch('data/rover.json')).json();
  for (const smp of ROVER.samples) {        // group samples cored at the same spot
    let site = SITES.find(x => x.lat === smp.lat && x.lon === smp.lon);
    if (!site) SITES.push(site = { lat: smp.lat, lon: smp.lon, samples: [] });
    site.samples.push(smp);
  }
  G.dLon = (G.b.east - G.b.west) / G.w; G.dLat = (G.b.north - G.b.south) / G.h;
}
const toIdx = ll => {
  const c = Math.floor((ll.lng - G.b.west) / G.dLon), r = Math.floor((G.b.north - ll.lat) / G.dLat);
  return (c < 0 || r < 0 || c >= G.w || r >= G.h) ? -1 : r * G.w + c;
};
const toLL = i => L.latLng(G.b.north - (Math.floor(i / G.w) + 0.5) * G.dLat, G.b.west + (i % G.w + 0.5) * G.dLon);

/* ---------- walking model ----------
 * Tobler's hiking function: v = 6·e^(−3.5·|S + 0.05|) km/h, S = rise/run.
 * Fastest on a gentle downhill, slower uphill and on steep descents. */
function moveSeconds(a, b, dist) {
  const S = (G.elev[b] - G.elev[a]) / dist;
  const kmh = 6 * Math.exp(-3.5 * Math.abs(S + 0.05)) * SUIT_FACTOR;
  return dist / (kmh / 3.6);
}
const NB = [[-1,-1],[0,-1],[1,-1],[-1,0],[1,0],[-1,1],[0,1],[1,1]];
const MAX_MS = 6 * SUIT_FACTOR / 3.6 * Math.exp(-3.5 * 0) ; // upper bound speed (m/s) for heuristic

/* ---------- binary min-heap on typed arrays ---------- */
class Heap {
  constructor(n) { this.k = new Float64Array(n); this.v = new Int32Array(n); this.n = 0; }
  push(key, val) { let i = this.n++; this.k[i] = key; this.v[i] = val;
    while (i > 0) { const p = (i - 1) >> 1; if (this.k[p] <= this.k[i]) break; this.swap(i, p); i = p; } }
  pop() { const val = this.v[0]; this.n--; this.k[0] = this.k[this.n]; this.v[0] = this.v[this.n];
    let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i;
      if (l < this.n && this.k[l] < this.k[m]) m = l; if (r < this.n && this.k[r] < this.k[m]) m = r;
      if (m === i) break; this.swap(i, m); i = m; } return val; }
  swap(a, b) { [this.k[a], this.k[b]] = [this.k[b], this.k[a]]; [this.v[a], this.v[b]] = [this.v[b], this.v[a]]; }
}

/* ---------- A* least-cost path ---------- */
function astar(start, goal, mult) {
  const N = G.w * G.h, g = new Float64Array(N).fill(Infinity), from = new Int32Array(N).fill(-1);
  const closed = new Uint8Array(N), heap = new Heap(N * 4);
  const gc = goal % G.w, gr = (goal / G.w) | 0;
  const hfun = i => Math.hypot(i % G.w - gc, ((i / G.w) | 0) - gr) * G.cell / MAX_MS;
  g[start] = 0; heap.push(hfun(start), start);
  while (heap.n) {
    const cur = heap.pop();
    if (cur === goal) break;
    if (closed[cur]) continue; closed[cur] = 1;
    const cc = cur % G.w, cr = (cur / G.w) | 0;
    for (const [dx, dy] of NB) {
      const nc = cc + dx, nr = cr + dy;
      if (nc < 0 || nr < 0 || nc >= G.w || nr >= G.h) continue;
      const nb = nr * G.w + nc, m = mult[G.haz[nb]];
      if (!isFinite(m) && nb !== goal) continue;
      const d = (dx && dy ? Math.SQRT2 : 1) * G.cell;
      const ng = g[cur] + moveSeconds(cur, nb, d) * (isFinite(m) ? m : 1);
      if (ng < g[nb]) { g[nb] = ng; from[nb] = cur; heap.push(ng + hfun(nb), nb); }
    }
  }
  if (from[goal] === -1 && goal !== start) return null;
  const path = []; for (let i = goal; i !== -1; i = from[i]) path.push(i);
  return path.reverse();
}

/* ---------- reverse Dijkstra: walking time from every cell back to the habitat ---------- */
function computeReturnTimes(home) {
  const N = G.w * G.h, t = new Float32Array(N).fill(Infinity), heap = new Heap(N * 4);
  const mult = PROFILES[0].mult; t[home] = 0; heap.push(0, home);
  const done = new Uint8Array(N);
  while (heap.n) {
    const u = heap.pop(); if (done[u]) continue; done[u] = 1;
    const uc = u % G.w, ur = (u / G.w) | 0;
    for (const [dx, dy] of NB) {
      const vc = uc + dx, vr = ur + dy;
      if (vc < 0 || vr < 0 || vc >= G.w || vr >= G.h) continue;
      const v = vr * G.w + vc; if (done[v] || G.haz[v] === NOGO) continue;
      const m = mult[G.haz[u]]; if (!isFinite(m) && u !== home) continue;
      const nt = t[u] + moveSeconds(v, u, (dx && dy ? Math.SQRT2 : 1) * G.cell) * m;  // walking v → u
      if (nt < t[v]) { t[v] = nt; heap.push(nt, v); }
    }
  }
  return t;
}

/* ---------- route statistics ---------- */
function routeStats(path) {
  let dist = 0, time = 0, gain = 0, loss = 0, maxSlope = 0, maxFromHome = 0;
  const hazM = [0, 0, 0, 0], cum = [{ d: 0, t: 0, z: G.elev[path[0]], h: G.haz[path[0]] }];
  const home = habitat.getLatLng();
  for (let k = 1; k < path.length; k++) {
    const a = path[k - 1], b = path[k];
    const d = ((Math.abs(a - b) === 1 || Math.abs(a - b) === G.w) ? 1 : Math.SQRT2) * G.cell;
    dist += d; time += moveSeconds(a, b, d);
    const dz = G.elev[b] - G.elev[a]; if (dz > 0) gain += dz; else loss -= dz;
    hazM[G.haz[b]] += d; maxSlope = Math.max(maxSlope, G.slope[b] / 10);
    maxFromHome = Math.max(maxFromHome, marsDist(home, toLL(b)));
    cum.push({ d: dist, t: time, z: G.elev[b], h: G.haz[b] });
  }
  // Turn-back point: last point where elapsed time + time to walk home ≤ EVA budget
  const budget = state.budgetH * 3600; let turnIdx = -1;
  for (let k = 0; k < path.length; k++) {
    if (cum[k].t + returnTime[path[k]] > budget) { turnIdx = Math.max(0, k - 1); break; }
  }
  const sites = SITES.filter(st => sitesNear(path, st));
  return { dist, time, gain, loss, maxSlope, maxFromHome, hazM, cum, turnIdx, sites,
           nSamples: sites.reduce((n, st) => n + st.samples.length, 0) };
}
// Distance on the Mars sphere (m)
function marsDist(a, b) {
  const R = 3396190, toR = Math.PI / 180;
  const x = (b.lng - a.lng) * toR * Math.cos((a.lat + b.lat) / 2 * toR), y = (b.lat - a.lat) * toR;
  return R * Math.hypot(x, y);
}

/* ---------- compute everything ---------- */
function recompute() {
  const home = toIdx(habitat.getLatLng());
  returnTime = computeReturnTimes(home);
  drawZone();
  routes = [];
  if (targets.length) {
    const base = [home, ...targets.map(t => toIdx(t.getLatLng())), home];
    for (const p of PROFILES) {
      let stops = base, r = routeThrough(stops, p);
      if (p.science && r.path) {
        // Greedy insertion: add the sample site with the smallest detour, keep it only
        // if the walk still fits the EVA budget (with 10% reserve) and the radius limit.
        const budget = state.budgetH * 3600 * 0.9, rad = state.radiusKm * 1000;
        const covered = path => SITES.filter(st => sitesNear(path, st));
        for (let added = 0; added < 4; added++) {
          const done = new Set(covered(r.path));
          let best = null;
          for (const st of SITES) {
            if (done.has(st)) continue;
            const ll = L.latLng(st.lat, st.lon), c = toIdx(ll);
            if (c < 0 || marsDist(habitat.getLatLng(), ll) > rad) continue;
            for (let k = 0; k < stops.length - 1; k++) {
              const A = toLL(stops[k]), B = toLL(stops[k + 1]);
              const extra = marsDist(A, ll) + marsDist(ll, B) - marsDist(A, B);
              if (!best || extra < best.extra) best = { extra, c, k };
            }
          }
          if (!best) break;
          const trial = [...stops.slice(0, best.k + 1), best.c, ...stops.slice(best.k + 1)];
          const tr = routeThrough(trial, p);
          if (!tr.path || tr.s.time > budget || tr.s.maxFromHome > rad) break;
          stops = trial; r = tr;
        }
        r.extraStops = stops.filter(i => !base.includes(i));
      }
      routes.push(r);
    }
  }
  drawRoutes(); renderCompare(); drawProfile(); renderTargets();
}

function routeThrough(stops, p) {
  let path = [];
  for (let k = 0; k < stops.length - 1; k++) {
    const leg = astar(stops[k], stops[k + 1], p.mult);
    if (!leg) return { p, path: null };
    path = path.concat(k ? leg.slice(1) : leg);
  }
  return { p, path, s: routeStats(path) };
}
// Is a sample site within 100 m of the route?
function sitesNear(path, st) {
  const ll = L.latLng(st.lat, st.lon);
  for (let k = 0; k < path.length; k += 2) if (marsDist(ll, toLL(path[k])) < 100) return true;
  return false;
}

/* ---------- drawing ---------- */
function drawZone() {
  if (zoneLayer) map.removeLayer(zoneLayer);
  if (!document.getElementById('showZone').checked) return;
  const cv = document.createElement('canvas'); cv.width = G.w; cv.height = G.h;
  const ctx = cv.getContext('2d'), img = ctx.createImageData(G.w, G.h);
  const half = state.budgetH * 3600 / 2, rad = state.radiusKm * 1000, home = habitat.getLatLng();
  for (let i = 0; i < G.w * G.h; i++) {
    const t = returnTime[i]; if (!(t <= half)) continue;
    if (marsDist(home, toLL(i)) > rad) continue;
    img.data[i * 4] = 62; img.data[i * 4 + 1] = 214; img.data[i * 4 + 2] = 224; img.data[i * 4 + 3] = 38;
  }
  ctx.putImageData(img, 0, 0);
  zoneLayer = L.imageOverlay(cv.toDataURL(), bounds(), { className: 'pixel', interactive: false }).addTo(map);
}
function drawRoutes() {
  routeLayers.forEach(l => map.removeLayer(l)); routeLayers = [];
  if (turnMarker) { map.removeLayer(turnMarker); turnMarker = null; }
  // radius ring drawn as polygon in degrees (circle in metres on Mars)
  const home = habitat.getLatLng(), ring = [];
  for (let a = 0; a <= 360; a += 4) {
    const r = state.radiusKm * 1000 / 3396190 * 180 / Math.PI;
    ring.push([home.lat + r * Math.sin(a * Math.PI / 180), home.lng + r * Math.cos(a * Math.PI / 180) / Math.cos(home.lat * Math.PI / 180)]);
  }
  routeLayers.push(L.polyline(ring, { color: '#3ed6e0', weight: 1, dashArray: '4 6', interactive: false }).addTo(map));
  for (const r of [...routes].sort((a, b) => (a.p.id === selected) - (b.p.id === selected))) {
    if (!r.path) continue;
    const sel = r.p.id === selected;
    const line = L.polyline(r.path.map(toLL), { color: r.p.color, weight: sel ? 4 : 2.5, opacity: sel ? 1 : 0.55,
      dashArray: sel ? null : '6 6' }).addTo(map);
    line.on('click', () => { selected = r.p.id; drawRoutes(); renderCompare(); drawProfile(); });
    routeLayers.push(line);
    if (r.extraStops) for (const i of r.extraStops) routeLayers.push(L.circleMarker(toLL(i),
      { radius: sel ? 9 : 6, color: r.p.color, weight: 2, fill: false, interactive: false }).addTo(map));
    if (sel && r.s.turnIdx >= 0) {
      turnMarker = L.marker(toLL(r.path[r.s.turnIdx]), { icon: pin('turn', '!') })
        .bindTooltip('Turn-back point: beyond here you cannot get home within the EVA budget').addTo(map);
    }
  }
}
const pin = (cls, txt) => L.divIcon({ className: '', html: `<div class="pin ${cls}">${txt}</div>`, iconSize: [26, 26], iconAnchor: [13, 13] });

/* ---------- panel ---------- */
const fmtT = s => `${Math.floor(s / 3600)}h ${String(Math.round(s % 3600 / 60)).padStart(2, '0')}m`;
function renderCompare() {
  const el = document.getElementById('compare');
  if (!routes.length) { el.innerHTML = '<p class="hint">Add a target to compute routes.</p>'; return; }
  const budget = state.budgetH * 3600, rad = state.radiusKm * 1000;
  const col = f => routes.map(r => `<td>${r.path ? f(r.s) : '<span class="bad">no route</span>'}</td>`).join('');
  const status = s => {
    if (s.time > budget) return '<span class="bad">OVER BUDGET</span>';
    if (s.maxFromHome > rad) return '<span class="bad">TOO FAR</span>';
    if (s.hazM[2] > 0) return '<span class="warn">CAUTION</span>';
    return '<span class="ok">GO</span>';
  };
  el.innerHTML = `<table><tr><th></th>${routes.map(r => `<th data-id="${r.p.id}" class="${r.p.id === selected ? 'sel' : ''}"><span class="dot" style="background:${r.p.color}"></span>${r.p.name}</th>`).join('')}</tr>
    <tr><td>Status</td>${col(status)}</tr>
    <tr><td>Sample sites passed</td>${col(s => s.sites.length ? `<b>${s.sites.length}</b> (${s.nSamples} cores)` : '0')}</tr>
    <tr><td>Distance</td>${col(s => (s.dist / 1000).toFixed(2) + ' km')}</tr>
    <tr><td>Walk time</td>${col(s => fmtT(s.time))}</tr>
    <tr><td>Climb ↑ / ↓</td>${col(s => `${Math.round(s.gain)} / ${Math.round(s.loss)} m`)}</tr>
    <tr><td>Max slope</td>${col(s => s.maxSlope.toFixed(1) + '°')}</tr>
    <tr><td>Moderate ground</td>${col(s => Math.round(s.hazM[1]) + ' m')}</tr>
    <tr><td>Difficult ground</td>${col(s => `<span class="${s.hazM[2] ? 'warn' : ''}">${Math.round(s.hazM[2])} m</span>`)}</tr>
    <tr><td>Farthest from home</td>${col(s => (s.maxFromHome / 1000).toFixed(2) + ' km')}</tr>
    <tr><td>Budget left</td>${col(s => `<span class="${s.time > budget ? 'bad' : 'ok'}">${s.time > budget ? '−' + fmtT(s.time - budget) : fmtT(budget - s.time)}</span>`)}</tr>
  </table>`;
  el.querySelectorAll('th[data-id]').forEach(th => th.onclick = () => { selected = th.dataset.id; drawRoutes(); renderCompare(); drawProfile(); });
}
function drawProfile() {
  const cv = document.getElementById('profile'), ctx = cv.getContext('2d');
  const W = cv.width = cv.clientWidth * devicePixelRatio, H = cv.height = 120 * devicePixelRatio;
  ctx.clearRect(0, 0, W, H);
  const r = routes.find(r => r.p.id === selected && r.path); if (!r) return;
  const c = r.s.cum, zs = c.map(p => p.z), zmin = Math.min(...zs), zmax = Math.max(...zs) + 1, D = c[c.length - 1].d;
  const pad = 22 * devicePixelRatio, X = d => pad + d / D * (W - pad * 1.4), Y = z => H - pad - (z - zmin) / (zmax - zmin) * (H - pad * 1.8);
  const hc = ['#59d98e', '#ffc83d', '#ff7a3d', '#ff4d5e'];
  ctx.lineWidth = 2 * devicePixelRatio;
  for (let k = 1; k < c.length; k++) { ctx.strokeStyle = hc[c[k].h]; ctx.beginPath(); ctx.moveTo(X(c[k - 1].d), Y(c[k - 1].z)); ctx.lineTo(X(c[k].d), Y(c[k].z)); ctx.stroke(); }
  if (r.s.turnIdx >= 0) { ctx.strokeStyle = '#ff4d5e'; ctx.setLineDash([4, 4]); const x = X(c[r.s.turnIdx].d); ctx.beginPath(); ctx.moveTo(x, 4); ctx.lineTo(x, H - pad); ctx.stroke(); ctx.setLineDash([]); }
  ctx.fillStyle = '#7c8496'; ctx.font = `${10 * devicePixelRatio}px JetBrains Mono`;
  ctx.fillText(`${Math.round(zmax)} m`, 4, 12 * devicePixelRatio); ctx.fillText(`${Math.round(zmin)} m`, 4, H - pad + 2);
  ctx.fillText('0 km', pad, H - 5); ctx.fillText(`${(D / 1000).toFixed(1)} km`, W - 50 * devicePixelRatio, H - 5);
  ctx.fillText(`Elevation profile · ${r.p.name}`, pad + 40 * devicePixelRatio, 12 * devicePixelRatio);
}
function renderTargets() {
  document.getElementById('targets').innerHTML = targets.map((t, i) => {
    const ll = t.getLatLng(); return `<li>T${i + 1} <span>${ll.lat.toFixed(4)}°N ${ll.lng.toFixed(4)}°E · ${Math.round(G.elev[toIdx(ll)])} m</span></li>`;
  }).join('');
}

function showCard(st) {
  const el = document.getElementById('card');
  el.innerHTML = st.samples.map(x => `<div class="sci">
      <div class="sci-h"><span class="num">#${x.num}</span> ${x.name}</div>
      <div class="sci-m">${x.rock} · rock “${x.feature}” · sealed sol ${x.sol}</div>
      <p>${x.why}</p></div>`).join('') +
    `<div class="sci-m">${st.lat.toFixed(4)}°N ${st.lon.toFixed(4)}°E · ${Math.round(G.elev[toIdx(L.latLng(st.lat, st.lon))])} m ·
     ${(marsDist(habitat.getLatLng(), L.latLng(st.lat, st.lon)) / 1000).toFixed(1)} km from habitat</div>
     <button id="addSite">✚ Add as target</button>`;
  document.getElementById('addSite').onclick = () => { addTarget(L.latLng(st.lat, st.lon)); recompute(); };
}

/* ---------- map setup ---------- */
const bounds = () => [[G.b.south, G.b.west], [G.b.north, G.b.east]];
function addTarget(ll) {
  if (toIdx(ll) < 0) return;
  const m = L.marker(ll, { icon: pin('tgt', 'T' + (targets.length + 1)), draggable: true }).addTo(map);
  m.on('dragend', recompute);
  m.on('contextmenu', () => { map.removeLayer(m); targets = targets.filter(t => t !== m); relabel(); recompute(); });
  targets.push(m);
}
function relabel() { targets.forEach((t, i) => t.setIcon(pin('tgt', 'T' + (i + 1)))); }

async function init() {
  await loadGrid();
  // Mars CRS: same lon/lat grid as EPSG:4326 but distances use Mars radius
  const MarsCRS = L.extend({}, L.CRS.EPSG4326, { R: 3396190 });
  map = L.map('map', { crs: MarsCRS, zoomSnap: 0.25, maxZoom: 18, attributionControl: true });
  map.fitBounds(bounds());
  const B = bounds();
  const layers = {
    'CTX image (MRO)': L.imageOverlay('data/ctx.jpg', B, { attribution: 'NASA/JPL/MSSS · USGS Astrogeology' }),
    'Elevation (CTX DTM)': L.imageOverlay('data/elevation.png', B),
    'Hillshade': L.imageOverlay('data/hillshade.png', B),
  };
  const overlays = {
    'Slope': L.imageOverlay('data/slope.png', B, { opacity: 0.8 }),
    'Hazard classes': L.imageOverlay('data/hazard.png', B, { className: 'pixel' }),
  };
  // Perseverance traverse + waypoints (Mars 2020 mission)
  const trav = L.layerGroup([L.polyline(ROVER.traverse.map(([x, y]) => [y, x]), { color: '#fff', weight: 1.5, opacity: 0.8, interactive: false })]);
  ROVER.waypoints.filter((_, i) => i % 6 === 0).forEach(([x, y, sol, km]) =>
    L.circleMarker([y, x], { radius: 2, color: '#fff', weight: 0, fillOpacity: 0.9 })
      .bindTooltip(`Perseverance · sol ${sol} · ${km} km driven`).addTo(trav));
  overlays[`Perseverance traverse (sol 0–${ROVER.last_sol})`] = trav;
  const sampleLayer = L.layerGroup();
  SITES.forEach(st => L.marker([st.lat, st.lon], { icon: L.divIcon({ className: '', html: '<div class="diamond"></div>', iconSize: [14, 14], iconAnchor: [7, 7] }) })
    .bindTooltip(st.samples.map(x => `#${x.num} ${x.name}`).join(' · '))
    .on('click', e => { L.DomEvent.stop(e); showCard(st); }).addTo(sampleLayer));
  overlays['Perseverance sample sites'] = sampleLayer;
  layers['CTX image (MRO)'].addTo(map); overlays['Hazard classes'].addTo(map); trav.addTo(map); sampleLayer.addTo(map);
  L.control.layers(layers, overlays, { collapsed: false }).addTo(map);
  L.control.scale({ imperial: false }).addTo(map);   // uses MarsCRS.R, so it's in Mars metres
  const lg = L.control({ position: 'bottomright' });
  lg.onAdd = () => { const d = L.DomUtil.create('div', 'legend'); d.innerHTML =
    '<i style="background:#ffc83d"></i>Moderate 8–15°<br><i style="background:#ff6e1e"></i>Difficult 15–25° / rough<br><i style="background:#e61e32"></i>No-go &gt;25°<br><i style="background:#3ed6e0;opacity:.5"></i>Safe-return zone'; return d; };
  lg.addTo(map);

  // Habitat at Perseverance landing site (Octavia E. Butler Landing)
  habitat = L.marker([18.4447, 77.4508], { icon: pin('hab', '⌂'), draggable: true }).addTo(map)
    .bindTooltip('Habitat (Octavia E. Butler Landing)');
  habitat.on('dragend', recompute);
  addTarget(L.latLng(18.4845, 77.4280));   // delta front (default demo target)

  map.on('click', e => {
    if (mode === 'habitat') { if (toIdx(e.latlng) >= 0) habitat.setLatLng(e.latlng); }
    else addTarget(e.latlng);
    recompute();
  });
  document.querySelectorAll('.mode').forEach(b => b.onclick = () => {
    mode = b.dataset.mode; document.querySelectorAll('.mode').forEach(x => x.classList.toggle('active', x === b));
    document.getElementById('hint').textContent = mode === 'habitat' ? 'Click the map to move the habitat.' :
      'Click the map to add science targets. Drag to move, right-click to delete.';
  });
  document.getElementById('clear').onclick = () => { targets.forEach(t => map.removeLayer(t)); targets = []; recompute(); };
  const bind = (id, key, unit) => { const el = document.getElementById(id), v = document.getElementById(id + 'V');
    const f = () => { state[key] = +el.value; v.textContent = el.value + unit; };
    f(); el.oninput = f; el.onchange = recompute; };
  bind('budget', 'budgetH', ' h'); bind('radius', 'radiusKm', ' km');
  document.getElementById('showZone').onchange = drawZone;
  document.getElementById('aboutLink').onclick = e => { e.preventDefault(); document.getElementById('about').classList.remove('hidden'); };
  document.getElementById('closeAbout').onclick = () => document.getElementById('about').classList.add('hidden');
  recompute();
}
init();
