// An interactive view of four human CFTR structures, drawn with raw WebGL2 (no library and no network call: the data is a module beside the page).
// The data are C-alpha atoms from RCSB PDB entries 6MSM, 5UAK, 8EIQ and 8EJ1 (CC0), superposed on 6MSM, and annotations from UniProt P13569 (CC BY 4.0);
// tools/cftr-model/ rebuilds them from the files. Every caption lives in the page's HTML so it can be read, reviewed and checked without this script; this file
// draws, and fills in numbers from the data.
//
// What is drawn, and what it does not claim:
//  - one smooth tube per run of neighbouring placed amino acids: a centripetal Catmull-Rom curve through the measured C-alpha positions, so it passes through
//    every measured point, but where it goes between two points is a drawing choice, not data. It is a real mesh with normals, lit and depth-tested;
//  - the tube widens into a flat band over the stretches UniProt annotates as a helix (wider) or a strand (narrower). That is UniProt's annotation, carried in the
//    data (uniprot.secondary) with the PDB ids UniProt cites; nothing here computes secondary structure from these coordinates. Turns and everything else are thin;
//  - colour by the part UniProt marks; darker where the chain is dense around it (a cheap stand-in for occlusion, computed from neighbour counts) and where it is
//    farther from the viewer. Shading is light on a shape, not a property of the protein;
//  - a dotted grey line across every stretch a model does not place (a straight connector, not the shape);
//  - residue 508 as a rose ring and dot, or an empty dashed ring midway between 507 and 509 where a structure lacks it;
//  - switching structure is a cross-dissolve of two complete pictures at a fixed camera. Nothing is interpolated between them, so nothing can look like
//    the protein moving from one state to another.
// Behaviour that matters for visitors: drawing happens only when something changed (or while the slow turn is on, at most 30 times a second), and never while the tab is
// hidden or the canvas is off the screen; the turn and the draw-in are off when the visitor asked for reduced motion, and the turn starts off on a page that sets
// data-turn="off" on #cftr-viewer; every drag has a slider or key alternative, and a lost WebGL context is rebuilt.
import DATA from './cftr-data.js';

const GREY = [0.55, 0.6, 0.66];
const ROSE = [1.0, 0.47, 0.59];
const WHITE = [1.0, 1.0, 1.0];
const PART_COLOURS = [[0.0, 0.95, 1.0], [0.82, 0.36, 1.0], [0.3, 0.55, 1.0], [1.0, 0.68, 0.1]];   // cyan, violet, blue, amber, in the order UniProt lists the four parts
const DOT_SPACING = 3.0;     // Angstroms between the dots that mark a stretch a model does not place
const FOV = 40 * Math.PI / 180;
const HOME = { theta: 0.6, phi: 1.35, zoom: 1 };
const HERO_ZOOM = 1.08;                   // the home page's stage is wider than tall; the protein is a little closer there
const TURN_RATE = 2 * Math.PI / 40;     // one slow revolution in forty seconds
const INTRO_MS = 2400;
const FADE_MS = 450;
const AMINO = { PHE: 'phenylalanine' };   // the data names residue 508 by its three-letter code; only the one residue the page talks about is spelled out

// ---- the tube: how finely the curve is cut and how thick each kind of stretch is drawn --------------------------------------------------------
const SUBDIV = 5;            // curve pieces between two neighbouring C-alpha atoms (3.8 Angstrom apart): about 0.76 Angstrom each
const SIDES = 8;             // vertices around the tube; phones draw this comfortably, and smooth normals hide the corners
// Cross-section half-sizes in Angstrom: [along the curvature normal, across it]. Equal values make a round tube; a helix is a wide flat band whose face follows
// its own coil, a strand a narrower band. The kind comes from UniProt's annotation per residue; between two residues the shape blends from one to the other.
const SHAPES = [[0.42, 0.42], [0.30, 1.25], [0.26, 0.95], [0.42, 0.42]];        // 0 everything else, 1 helix, 2 strand, 3 turn
const SEC_CODE = { helix: 1, strand: 2, turn: 3 };
const DENSE_RADIUS = 9.0;    // Angstrom: neighbours within this distance count towards the darkening of crowded stretches

const VS_TUBE = `#version 300 es
in vec3 aPos; in vec3 aNrm; in float aRes; in float aAo;
uniform mat4 uMvp; uniform mat4 uView;
out vec3 vN; out vec3 vP; out float vRes; out float vAo;
void main(){
  gl_Position = uMvp * vec4(aPos, 1.0);
  vP = (uView * vec4(aPos, 1.0)).xyz;
  vN = mat3(uView) * aNrm;
  vRes = aRes; vAo = aAo;
}`;
const FS_TUBE = `#version 300 es
precision mediump float;
in vec3 vN; in vec3 vP; in float vRes; in float vAo;
uniform sampler2D uPal; uniform int uPalMax;
uniform float uAlpha; uniform float uGrey; uniform vec3 uGreyColour;
uniform float uProgress; uniform float uWalk; uniform float uWalkK; uniform float uComet;
uniform float uDist; uniform float uRad;
out vec4 o;
const vec3 KEY = vec3(-0.4558, 0.7090, 0.5381);     // a soft light from the upper left, a little in front of the viewer (unit length)
const vec3 FILL = vec3(0.8138, -0.3488, 0.4650);    // a dim fill from the lower right
void main(){
  if (vRes > uProgress) discard;                                             // the draw-in: only the chain up to this residue exists yet
  vec3 n = normalize(vN);
  if (!gl_FrontFacing) n = -n;
  vec3 v = normalize(-vP);
  vec3 base = texelFetch(uPal, ivec2(clamp(int(floor(vRes + 0.5)), 0, uPalMax), 0), 0).rgb;
  base = mix(base, uGreyColour, uGrey);
  float occ = 1.0 - 0.6 * vAo;                                               // crowded stretches sit in shadow
  float kd = dot(n, KEY) * 0.5 + 0.5; kd = kd * kd * kd;                     // wrap-around key light with a soft fall-off: no hard terminator on a thin tube
  float fd = max(dot(n, FILL), 0.0);
  float spec = pow(max(dot(n, normalize(KEY + v)), 0.0), 48.0) * 0.45 * occ;
  float rim = pow(1.0 - max(dot(n, v), 0.0), 3.5);                           // a neon edge where the surface turns away
  vec3 col = base * (0.1 + 0.9 * kd) * occ + base * 0.22 * fd * occ + spec * mix(vec3(1.0), base, 0.25) + rim * (base * 0.8 + vec3(0.05, 0.04, 0.1));
  float near = clamp(0.5 + (uDist + vP.z) / (2.0 * uRad) * 0.9, 0.0, 1.0);  // vP.z is negative in front of the camera
  col = mix(col, vec3(0.02, 0.02, 0.045), 0.5 * (1.0 - near));              // the far side recedes: darker, not thinner or dimmer in any data sense
  if (uWalk > 0.0) {
    float t = uWalk - vRes, f = 1.0;
    if (t < 0.0) f = 0.7;
    else if (t < uComet) f = 1.0 + 0.9 * (1.0 - t / uComet);
    else f = 0.85;
    col *= mix(1.0, f, uWalkK);
  }
  o = vec4(min(col, vec3(1.0)) * uAlpha, uAlpha);
}`;
const VS_POINT = `#version 300 es
in vec3 aPos;
uniform mat4 uMvp; uniform vec3 uPos; uniform float uSize; uniform float uUseAttr;
void main(){ gl_Position = uMvp * vec4(uUseAttr > 0.5 ? aPos : uPos, 1.0); gl_PointSize = uSize; }`;
const FS_POINT = `#version 300 es
precision mediump float;
uniform vec3 uColour; uniform float uMode; uniform float uAlpha;
out vec4 o;
void main(){
  float d = length(gl_PointCoord - 0.5);
  if (d > 0.5) discard;
  float a;
  if (uMode > 2.5) {                                            // dashed ring: nothing is placed here
    if (d < 0.36) discard;
    float ang = atan(gl_PointCoord.y - 0.5, gl_PointCoord.x - 0.5) / 6.2831853 + 0.5;
    if (fract(ang * 10.0) > 0.55) discard;
    a = 0.95;
  }
  else if (uMode > 1.5) a = smoothstep(0.5, 0.0, d) * 0.55;    // halo
  else if (uMode > 0.5) { if (d < 0.36) discard; a = 0.95; }    // ring
  else a = smoothstep(0.5, 0.2, d);                              // dot
  a *= uAlpha;
  o = vec4(uColour * a, a);
}`;

// ---- small matrix helpers (column-major, as WebGL wants) ------------------------------------------------------------------------------
function perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2), nf = 1 / (near - far);
  return [f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0];
}
function lookAt(eye, up) {                   // a camera at `eye` looking at the origin; `up` must not be parallel to the line of sight
  const unit = (v) => { const l = Math.hypot(v[0], v[1], v[2]); return [v[0] / l, v[1] / l, v[2] / l]; };
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const z = unit(eye), x = unit(cross(up, z)), y = cross(z, x);
  return [x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0, -dot(x, eye), -dot(y, eye), -dot(z, eye), 1];
}
function multiply(a, b) {
  const o = new Array(16).fill(0);
  for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) for (let k = 0; k < 4; k++) o[c * 4 + r] += a[k * 4 + r] * b[c * 4 + k];
  return o;
}

// ---- reading the data ---------------------------------------------------------------------------------------------------------------
function partsOf(data) {
  const kinds = { membrane: 0, binding: 0 };
  return data.domains.filter(d => d.type === 'Domain').map((d, i) => {
    const membrane = /transmembrane/i.test(d.name);
    const n = ++kinds[membrane ? 'membrane' : 'binding'];
    return { start: d.start, end: d.end, colour: PART_COLOURS[i % PART_COLOURS.length], label: (membrane ? 'Membrane-spanning part ' : 'Nucleotide-binding part ') + n };
  });
}
function colourOf(parts, n) {
  const p = parts.find(q => n >= q.start && n <= q.end);
  return p ? p.colour : GREY;
}
function indexStructure(st, length) {
  const covered = new Uint8Array(length + 2), at = new Map();
  for (const seg of st.segments) {
    for (let i = 0; i < seg.xyz.length / 3; i++) {
      covered[seg.start + i] = 1;
      at.set(seg.start + i, [seg.xyz[3 * i], seg.xyz[3 * i + 1], seg.xyz[3 * i + 2]]);
    }
  }
  st.covered = covered; st.at = at;
  return st;
}
function secondaryByResidue(data) {        // UniProt's helix / strand / turn annotation as one code per residue number (0 where it marks nothing)
  const t = new Uint8Array(data.length + 2);
  for (const s of data.uniprot.secondary || []) for (let n = s.start; n <= s.end && n < t.length; n++) t[n] = SEC_CODE[s.type] || 0;
  return t;
}

// ---- the mesh: built once per structure on the CPU, kept so a lost context only needs a re-upload ------------------------------------------
function crowding(st, length) {
  // How many other placed C-alpha atoms lie within DENSE_RADIUS of each one, scaled to 0..1. A cheap stand-in for ambient occlusion: a lobe's core is crowded, a loop on the
  // surface is not. Sorted by x so each pair is tested once and most are skipped early.
  const pts = [];
  for (const seg of st.segments) for (let i = 0; i < seg.xyz.length / 3; i++) pts.push([seg.start + i, seg.xyz[3 * i], seg.xyz[3 * i + 1], seg.xyz[3 * i + 2]]);
  pts.sort((a, b) => a[1] - b[1]);
  const count = new Float32Array(length + 2), r2 = DENSE_RADIUS * DENSE_RADIUS;
  for (let i = 0; i < pts.length; i++) {
    const p = pts[i];
    for (let j = i + 1; j < pts.length; j++) {
      const q = pts[j], dx = q[1] - p[1];
      if (dx > DENSE_RADIUS) break;
      const dy = q[2] - p[2], dz = q[3] - p[3];
      if (dx * dx + dy * dy + dz * dz < r2) { count[p[0]]++; count[q[0]]++; }
    }
  }
  const ao = new Float32Array(length + 2);
  for (let n = 0; n < ao.length; n++) ao[n] = Math.max(0, Math.min(1, (count[n] - 10) / 20));
  return ao;
}

function tubeMesh(st, sec, ao) {
  // One tube per run of consecutive placed residues. Rings of SIDES vertices sit on a centripetal Catmull-Rom curve through the C-alpha positions, SUBDIV rings per
  // residue step, plus a flat cap at each end of a run. Vertex: position (3), normal (3), residue number as a float (1), crowding (1): 8 floats.
  let rings = 0, runs = 0;
  for (const seg of st.segments) { const n = seg.xyz.length / 3; if (n < 2) continue; rings += (n - 1) * SUBDIV + 1; runs++; }
  const nVerts = rings * SIDES + runs * 2 * (SIDES + 1);
  const nTris = (rings - runs) * SIDES * 2 + runs * 2 * SIDES;
  const verts = new Float32Array(nVerts * 8);
  const idx = nVerts > 65535 ? new Uint32Array(nTris * 3) : new Uint16Array(nTris * 3);
  const cosT = [], sinT = [];
  for (let j = 0; j < SIDES; j++) { cosT.push(Math.cos(2 * Math.PI * j / SIDES)); sinT.push(Math.sin(2 * Math.PI * j / SIDES)); }
  let v = 0, f = 0;
  const put = (x, y, z, nx, ny, nz, res, occ) => { verts[v++] = x; verts[v++] = y; verts[v++] = z; verts[v++] = nx; verts[v++] = ny; verts[v++] = nz; verts[v++] = res; verts[v++] = occ; };
  const tri = (a, b, c) => { idx[f++] = a; idx[f++] = b; idx[f++] = c; };
  const T = [0, 0, 0], N = [0, 0, 0], B = [0, 0, 0], prevN = [0, 0, 1];
  for (const seg of st.segments) {
    const n = seg.xyz.length / 3;
    if (n < 2) continue;
    const P = seg.xyz;
    // one direction hint per residue: the discrete second difference, which points at a helix's own axis and alternates across a strand's pleats
    const H = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) {
      const a = Math.max(1, Math.min(n - 2, i));                                   // the first and last residue borrow their neighbour's hint
      let hx = P[3 * (a - 1)] - 2 * P[3 * a] + P[3 * (a + 1)], hy = P[3 * (a - 1) + 1] - 2 * P[3 * a + 1] + P[3 * (a + 1) + 1], hz = P[3 * (a - 1) + 2] - 2 * P[3 * a + 2] + P[3 * (a + 1) + 2];
      const l = Math.hypot(hx, hy, hz);
      if (l < 1e-4) { if (i > 0) { hx = H[3 * i - 3]; hy = H[3 * i - 2]; hz = H[3 * i - 1]; } else { hx = 0; hy = 0; hz = 1; } }
      else { hx /= l; hy /= l; hz /= l; }
      // keep the hint on one side from residue to residue, except inside a helix, where it is meant to turn around the axis (about 100 degrees a residue)
      if (i > 0 && sec[seg.start + i] !== 1 && hx * H[3 * i - 3] + hy * H[3 * i - 2] + hz * H[3 * i - 1] < 0) { hx = -hx; hy = -hy; hz = -hz; }
      H[3 * i] = hx; H[3 * i + 1] = hy; H[3 * i + 2] = hz;
    }
    // the curve: ring positions, then tangents by central difference
    const pos = splinePositions(P);
    const R = pos.length / 3;
    const base = v / 8;
    for (let r = 0; r < R; r++) {
      const i = Math.min(Math.floor(r / SUBDIV), n - 2), u = (r - i * SUBDIV) / SUBDIV;
      const ra = Math.max(0, r - 1), rb = Math.min(R - 1, r + 1);
      T[0] = pos[3 * rb] - pos[3 * ra]; T[1] = pos[3 * rb + 1] - pos[3 * ra + 1]; T[2] = pos[3 * rb + 2] - pos[3 * ra + 2];
      normalise(T);
      N[0] = H[3 * i] * (1 - u) + H[3 * i + 3] * u; N[1] = H[3 * i + 1] * (1 - u) + H[3 * i + 4] * u; N[2] = H[3 * i + 2] * (1 - u) + H[3 * i + 5] * u;
      if (!perpendicular(N, T)) { N[0] = prevN[0]; N[1] = prevN[1]; N[2] = prevN[2]; if (!perpendicular(N, T)) anyPerpendicular(N, T); }
      prevN[0] = N[0]; prevN[1] = N[1]; prevN[2] = N[2];
      B[0] = T[1] * N[2] - T[2] * N[1]; B[1] = T[2] * N[0] - T[0] * N[2]; B[2] = T[0] * N[1] - T[1] * N[0];
      const s0 = SHAPES[sec[seg.start + i]], s1 = SHAPES[sec[seg.start + i + 1]];
      const a = s0[0] * (1 - u) + s1[0] * u, b = s0[1] * (1 - u) + s1[1] * u;
      const res = seg.start + i + u, occ = ao[seg.start + i] * (1 - u) + ao[seg.start + i + 1] * u;
      const cx = pos[3 * r], cy = pos[3 * r + 1], cz = pos[3 * r + 2];
      if (r === 0) {                                                             // cap at the start of the run, facing back along the curve
        const c = base + r * SIDES;
        put(cx, cy, cz, -T[0], -T[1], -T[2], res, occ);
        for (let j = 0; j < SIDES; j++) put(cx + N[0] * a * cosT[j] + B[0] * b * sinT[j], cy + N[1] * a * cosT[j] + B[1] * b * sinT[j], cz + N[2] * a * cosT[j] + B[2] * b * sinT[j], -T[0], -T[1], -T[2], res, occ);
        for (let j = 0; j < SIDES; j++) tri(c, c + 1 + (j + 1) % SIDES, c + 1 + j);
      }
      const ring = v / 8;
      for (let j = 0; j < SIDES; j++) {
        const c = cosT[j], s = sinT[j];
        const nx = N[0] * c / a + B[0] * s / b, ny = N[1] * c / a + B[1] * s / b, nz = N[2] * c / a + B[2] * s / b, l = Math.hypot(nx, ny, nz);
        put(cx + N[0] * a * c + B[0] * b * s, cy + N[1] * a * c + B[1] * b * s, cz + N[2] * a * c + B[2] * b * s, nx / l, ny / l, nz / l, res, occ);
      }
      if (r > 0) {
        const prev = ring - SIDES;
        for (let j = 0; j < SIDES; j++) {
          const j1 = (j + 1) % SIDES;
          tri(prev + j, ring + j, prev + j1);
          tri(prev + j1, ring + j, ring + j1);
        }
      }
      if (r === R - 1) {                                                         // cap at the end of the run, facing on along the curve
        const c = v / 8;
        put(cx, cy, cz, T[0], T[1], T[2], res, occ);
        for (let j = 0; j < SIDES; j++) put(cx + N[0] * a * cosT[j] + B[0] * b * sinT[j], cy + N[1] * a * cosT[j] + B[1] * b * sinT[j], cz + N[2] * a * cosT[j] + B[2] * b * sinT[j], T[0], T[1], T[2], res, occ);
        for (let j = 0; j < SIDES; j++) tri(c, c + 1 + j, c + 1 + (j + 1) % SIDES);
      }
    }
  }
  return { verts: verts.subarray(0, v), idx: idx.subarray(0, f), triangles: f / 3, vertices: v / 8 };
}
function splinePositions(P) {
  // Ring centres along a centripetal Catmull-Rom curve through one run of measured C-alpha positions (P is the run as a flat x,y,z list, two or more points).
  // SUBDIV rings sit between neighbouring points, and every SUBDIV-th ring is exactly a measured point: the curve passes through all of them.
  const n = P.length / 3;
  const px = (i, k) => (i < 0 ? 2 * P[k] - P[3 + k] : i >= n ? 2 * P[3 * (n - 1) + k] - P[3 * (n - 2) + k] : P[3 * i + k]);    // ends extend straight on
  const R = (n - 1) * SUBDIV + 1;
  const pos = new Float32Array(R * 3);
  const knot = (a, b) => Math.sqrt(Math.hypot(px(a, 0) - px(b, 0), px(a, 1) - px(b, 1), px(a, 2) - px(b, 2))) || 1e-3;
  for (let r = 0; r < R; r++) {
    const i = Math.min(Math.floor(r / SUBDIV), n - 2), u = (r - i * SUBDIV) / SUBDIV;
    const t1 = knot(i - 1, i), t2 = t1 + knot(i, i + 1), t3 = t2 + knot(i + 1, i + 2);
    for (let k = 0; k < 3; k++) pos[3 * r + k] = catmullRom(px(i - 1, k), px(i, k), px(i + 1, k), px(i + 2, k), t1, t2, t3, t1 + (t2 - t1) * u);
  }
  return pos;
}
function catmullRom(p0, p1, p2, p3, t1, t2, t3, t) {
  // One coordinate of a centripetal Catmull-Rom curve (Barry and Goldman's form) between p1 (knot t1) and p2 (knot t2); the first knot is 0. Knot spacing is the
  // square root of the distance between control points, which keeps the curve from looping or overshooting between close points.
  const A1 = ((t1 - t) * p0 + t * p1) / t1, A2 = ((t2 - t) * p1 + (t - t1) * p2) / (t2 - t1), A3 = ((t3 - t) * p2 + (t - t2) * p3) / (t3 - t2);
  const B1 = ((t2 - t) * A1 + t * A2) / t2, B2 = ((t3 - t) * A2 + (t - t1) * A3) / (t3 - t1);
  return ((t2 - t) * B1 + (t - t1) * B2) / (t2 - t1);
}
function normalise(a) { const l = Math.hypot(a[0], a[1], a[2]) || 1; a[0] /= l; a[1] /= l; a[2] /= l; return a; }
function perpendicular(nrm, t) {            // make `nrm` perpendicular to the unit vector `t`, in place; false when nothing of it is left
  const d = nrm[0] * t[0] + nrm[1] * t[1] + nrm[2] * t[2];
  nrm[0] -= d * t[0]; nrm[1] -= d * t[1]; nrm[2] -= d * t[2];
  const l = Math.hypot(nrm[0], nrm[1], nrm[2]);
  if (l < 1e-3) return false;
  nrm[0] /= l; nrm[1] /= l; nrm[2] /= l;
  return true;
}
function anyPerpendicular(nrm, t) {
  const ax = Math.abs(t[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  nrm[0] = ax[0]; nrm[1] = ax[1]; nrm[2] = ax[2];
  perpendicular(nrm, t);
}

const meshCache = new Map();               // structure index -> { verts, idx, triangles, vertices }: the CPU side survives a lost context
function meshes(data, sec) {
  const t0 = performance.now();
  if (!meshCache.size) data.structures.forEach((st, i) => meshCache.set(i, tubeMesh(st, sec, crowding(st, data.length))));
  return { built: meshCache, ms: performance.now() - t0 };
}

function build(gl, data, colourFn, sec) {
  const compile = (type, src) => {
    const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error('shader: ' + gl.getShaderInfoLog(s));
    return s;
  };
  const program = (vs, fs, attrs) => {
    const p = gl.createProgram();
    gl.attachShader(p, compile(gl.VERTEX_SHADER, vs)); gl.attachShader(p, compile(gl.FRAGMENT_SHADER, fs));
    attrs.forEach((a, i) => gl.bindAttribLocation(p, i, a));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error('link: ' + gl.getProgramInfoLog(p));
    return p;
  };
  const uniforms = (p, names) => Object.fromEntries(names.map(k => [k, gl.getUniformLocation(p, k)]));
  const t0 = performance.now();
  const { built, ms: meshMs } = meshes(data, sec);
  const gpu = { tube: program(VS_TUBE, FS_TUBE, ['aPos', 'aNrm', 'aRes', 'aAo']), point: program(VS_POINT, FS_POINT, ['aPos']), structures: [], meshMs, palette: null, triangles: [] };
  gpu.tubeU = uniforms(gpu.tube, ['uMvp', 'uView', 'uPal', 'uPalMax', 'uAlpha', 'uGrey', 'uGreyColour', 'uProgress', 'uWalk', 'uWalkK', 'uComet', 'uDist', 'uRad']);
  gpu.pointU = uniforms(gpu.point, ['uMvp', 'uPos', 'uSize', 'uUseAttr', 'uColour', 'uMode', 'uAlpha']);
  data.structures.forEach((st, i) => {
    const mesh = built.get(i);
    const vao = gl.createVertexArray(), buf = gl.createBuffer(), ibuf = gl.createBuffer();
    gl.bindVertexArray(vao);
    gl.bindBuffer(gl.ARRAY_BUFFER, buf); gl.bufferData(gl.ARRAY_BUFFER, mesh.verts, gl.STATIC_DRAW);
    for (const [loc, size, off] of [[0, 3, 0], [1, 3, 12], [2, 1, 24], [3, 1, 28]]) { gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 32, off); }
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ibuf); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, mesh.idx, gl.STATIC_DRAW);
    gl.bindVertexArray(null);
    const dots = [];                          // a dotted straight connector across each stretch the model does not place: a hint of the gap, not its shape
    for (let s = 0; s + 1 < st.segments.length; s++) {
      const a = st.segments[s].xyz.slice(-3), b = st.segments[s + 1].xyz.slice(0, 3);
      const len = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
      const k = Math.max(1, Math.floor(len / DOT_SPACING));
      for (let j = 1; j <= k; j++) {
        const t = j / (k + 1);
        dots.push(a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t);
      }
    }
    const dvao = gl.createVertexArray(), dbuf = gl.createBuffer();
    gl.bindVertexArray(dvao); gl.bindBuffer(gl.ARRAY_BUFFER, dbuf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(dots), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 12, 0);
    gl.bindVertexArray(null);
    gpu.structures.push({ vao, count: mesh.idx.length, type: mesh.idx instanceof Uint32Array ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT, dots: dvao, dotCount: dots.length / 3 });
    gpu.triangles.push(mesh.triangles);
  });
  // colour by residue number lives in a one-row texture, so recolouring (a topology layer, say) rewrites 1,481 texels and no geometry
  gpu.palette = gl.createTexture(); gpu.paletteN = data.length + 1;
  gl.bindTexture(gl.TEXTURE_2D, gpu.palette);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, gpu.paletteN, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
  recolour(gl, gpu, colourFn);
  gpu.buildMs = performance.now() - t0;
  return gpu;
}

function recolour(gl, gpu, colourFn) {
  const px = new Uint8Array(gpu.paletteN * 4);
  for (let n = 1; n < gpu.paletteN; n++) {
    const c = colourFn(n);
    px[4 * n] = Math.round(Math.max(0, Math.min(1, c[0])) * 255); px[4 * n + 1] = Math.round(Math.max(0, Math.min(1, c[1])) * 255); px[4 * n + 2] = Math.round(Math.max(0, Math.min(1, c[2])) * 255); px[4 * n + 3] = 255;
  }
  gl.bindTexture(gl.TEXTURE_2D, gpu.palette);
  gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, gpu.paletteN, 1, gl.RGBA, gl.UNSIGNED_BYTE, px);
}

// ---- the viewer ---------------------------------------------------------------------------------------------------------------------
export async function start() {
  const listeners = [];
  const layerFns = new Map();
  const canvas = document.getElementById('cftr-canvas');
  if (!canvas) return;
  const fallback = document.getElementById('cftr-fallback');
  const root = document.getElementById('cftr-viewer');
  const hero = !!root && root.dataset.mode === 'hero';
  const fail = (why) => {
    if (fallback) { fallback.hidden = false; const w = fallback.querySelector('[data-why]'); if (w) w.textContent = why; }
    if (root) root.hidden = true;
  };
  const data = DATA;
  const gl = canvas.getContext('webgl2', { antialias: true, alpha: true, premultipliedAlpha: true, depth: true });
  if (!gl) return fail('Your browser could not start WebGL2.');

  const parts = partsOf(data);
  const length = data.length;
  const sec = secondaryByResidue(data);
  data.structures.forEach(st => indexStructure(st, length));
  let radius = 1;
  for (const st of data.structures) for (const p of st.at.values()) radius = Math.max(radius, Math.hypot(p[0], p[1], p[2]));
  const byId = Object.fromEntries(data.structures.map((s, i) => [s.id, i]));
  const reduceQuery = window.matchMedia('(prefers-reduced-motion: reduce)');
  const turnOff = !!root && root.dataset.turn === 'off';                   // a page can ask for the slow turn to start off (the home page does)

  const home = hero ? { ...HOME, zoom: HERO_ZOOM } : HOME;
  const S = { i: 0, overlay: false, walk: null, theta: home.theta, phi: home.phi, zoom: home.zoom, playing: !reduceQuery.matches && !turnOff, fade: null, intro: null, layers: new Set(), viewTouched: false, walkAt: 0 };
  let colourFn = (n) => colourOf(parts, n);
  let hashTimer = 0;
  const draws = { tube: 0, dots: 0, depthOn: null, absentRingMode: null };       // what was drawn, published on the canvas for tests
  let gpu = null, raf = 0, lastT = 0, lastDrawT = 0, dirty = true, lost = false, offscreen = false, frames = 0, w = 1, h = 1, dpr = 1, seenW = 0, seenH = 0;
  const $ = (id) => document.getElementById(id);
  const el = { readout: $('cftr-readout'), walkOut: $('cftr-walk-out'), walk: $('cftr-walk'), turn: $('cftr-turn'), tilt: $('cftr-tilt'), play: $('cftr-play'),
               overlay: $('cftr-overlay'), legend: $('cftr-legend'), picks: $('cftr-picks'), caption: $('cftr-caption'), note: $('cftr-note'), ruler: $('cftr-ruler-svg') };
  const HASH_MODE = root ? (root.dataset.hash || '') : '';
  const fullLink = document.getElementById('cftr-full-link');
  const on = (node, type, fn, opts) => { if (node) node.addEventListener(type, fn, opts); };

  const residueName = (st, n) => (st.res508 && n === 508) ? (AMINO[st.res508.name] || st.res508.name.toLowerCase()) : null;
  const partName = (n) => { const p = parts.find(q => n >= q.start && n <= q.end); return p ? p.label.toLowerCase() : (n >= 654 && n <= 831 ? 'the stretch that UniProt marks as disordered (654 to 831)' : 'no marked part'); };
  const stretchAround = (st, n) => {
    let a = n, b = n;
    while (a > 1 && !st.covered[a - 1]) a--;
    while (b < length && !st.covered[b + 1]) b++;
    return [a, b];
  };

  // ---- text that depends on the data ----
  function walkText() {
    const st = data.structures[S.i];
    if (S.walk === null) return hero ? '' : 'No residue is marked. Move the slider to walk along the chain.';
    const n = S.walk;
    if (st.covered[n]) {
      const nm = residueName(st, n);
      return 'Residue ' + n + (nm ? ' (' + nm + ')' : '') + ': placed in this model, in ' + partName(n) + '.';
    }
    if (n === 508 && st.covered[507] && st.covered[509]) return 'Residue 508 is absent from this structure: the model runs from residue 507 straight to 509.';
    const [a, b] = stretchAround(st, n);
    return 'Residue ' + n + ': not placed in this model (the stretch ' + a + ' to ' + b + ' is not modelled), in ' + partName(n) + '.';
  }
  function captionText(st) {
    const r508 = st.res508 ? 'Residue 508 is present in this structure.' : 'Residue 508 is absent from this structure.';
    return 'CFTR is a human protein, and variants of it are linked to cystic fibrosis. Structure ' + st.id + ' places ' + st.residues_traced.toLocaleString('en-US') + ' of its '
      + length.toLocaleString('en-US') + ' amino acids. ' + r508 + ' A model fitted to electron-microscope data, not a photograph.';
  }
  function updateReadout() {
    const st = data.structures[S.i];
    if (el.readout) {
      const lines = [
        [st.id, st.title],
        ['Method', 'electron microscopy, ' + st.resolution_A + ' Å; deposited ' + st.deposited],
        ['Model', st.residues_traced + ' of ' + length + ' residues are placed.'],
        ['Residue 508', st.res508 ? 'present in this structure (' + residueName(st, 508) + ')' : 'absent from this structure: the chain runs from 507 to 509'],
      ];
      if (S.overlay && S.i !== 0) lines.push(['Lined up on 6MSM', st.superposition_residues + ' membrane-spanning residues, ' + st.superposition_rmsd_A.toFixed(2) + ' Å apart (computed for this page)']);
      el.readout.replaceChildren(...lines.map(([k, v]) => { const p = document.createElement('p'); const b = document.createElement('strong'); b.textContent = k + ': '; p.append(b, v); return p; }));
    }
    if (el.caption) el.caption.textContent = captionText(st);
    if (el.walkOut) el.walkOut.textContent = walkText();
    if (el.walk) el.walk.setAttribute('aria-valuetext', S.walk === null ? 'Residue ' + el.walk.value + ', no marker' : walkText());
    drawRuler();
    drawSequence();
    publish();
  }
  function publish() {                       // the current state as data on the canvas, for tests and for anyone inspecting the page
    canvas.dataset.state = JSON.stringify({ structure: data.structures[S.i].id, overlay: S.overlay, walk: S.walk, playing: S.playing, fading: !!S.fade, intro: !!S.intro, offscreen, layers: [...S.layers] });
    for (const fn of listeners) { try { fn(snapshot()); } catch (e) { /* a listener must not stop the viewer */ } }
    writeHashSoon();
  }
  function snapshot() {
    return { structure: data.structures[S.i].id, overlay: S.overlay, walk: S.walk, theta: S.theta, phi: S.phi, zoom: S.zoom, playing: S.playing, layers: [...S.layers], viewTouched: S.viewTouched };
  }
  // ---- a view as a short address: structure, marker, overlay, camera, layers. Nothing leaves the browser; it is only a link the visitor can copy.
  function encodeHash() {
    const q = ['s=' + data.structures[S.i].id];
    if (S.walk !== null) q.push('r=' + S.walk);
    if (S.overlay && S.i !== 0) q.push('o=1');
    if (S.viewTouched) q.push('y=' + Math.round(((S.theta * 180 / Math.PI) % 360 + 360) % 360), 'p=' + Math.round(S.phi * 180 / Math.PI), 'z=' + S.zoom.toFixed(2));
    if (S.layers.size) q.push('l=' + [...S.layers].sort().join(','));
    return q.join('&');
  }
  function writeHashSoon() {
    if (!root || !HASH_MODE) return;
    clearTimeout(hashTimer);
    hashTimer = setTimeout(() => {
      const h = encodeHash();
      if (HASH_MODE === 'url') { try { history.replaceState(null, '', '#' + h); } catch (e) { /* not allowed here */ } }
      else if (fullLink) fullLink.setAttribute('href', 'cf/cftr/#' + h);
    }, 200);
  }

  const sameColour = (a, b) => a === b || (a[0] === b[0] && a[1] === b[1] && a[2] === b[2]);
  // ---- UniProt's reference sequence around the marked residue: lit where this structure places it, an empty box where the file says it was deleted
  const AMINO_ONE = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val' };
  function drawSequence() {
    const box = $('cftr-seq');
    if (!box) return;
    const st = data.structures[S.i], centre = S.walk ?? 508, seq = data.uniprot.sequence;
    const deleted = new Set(st.differences.filter(d => d.details !== 'engineered mutation' && d.position && !st.covered[d.position] && !d.file_residue).map(d => d.position));
    const spans = [];
    for (let n = Math.max(1, centre - 10); n <= Math.min(length, centre + 10); n++) {
      const s = document.createElement('span');
      s.textContent = deleted.has(n) ? ' ' : seq[n - 1];
      s.className = 'cftr-seq-cell ' + (deleted.has(n) ? 'is-absent' : st.covered[n] ? 'is-placed' : 'is-unplaced') + (n === centre ? ' is-centre' : '');
      s.title = deleted.has(n) ? n + ': the file says this residue is deleted' : n + ' ' + AMINO_ONE[seq[n - 1]] + (st.covered[n] ? '' : ' (not placed in this structure)');
      spans.push(s);
    }
    box.replaceChildren(...spans);
  }
  // ---- the chain ruler: the whole chain as a strip, coloured like the chain, with the places a model leaves out ----
  const SVGNS = 'http://www.w3.org/2000/svg';
  function drawRuler() {
    if (!el.ruler) return;
    const st = data.structures[S.i], kids = [];
    const rect = (x, wd, y, ht, fill, extra) => {
      const r = document.createElementNS(SVGNS, 'rect');
      r.setAttribute('x', String(x)); r.setAttribute('width', String(wd)); r.setAttribute('y', String(y)); r.setAttribute('height', String(ht)); r.setAttribute('fill', fill);
      if (extra) for (const [k, v] of Object.entries(extra)) r.setAttribute(k, v);
      return r;
    };
    const css = (c) => 'rgb(' + c.map(v => Math.round(v * 255)).join(',') + ')';
    let n = 1;
    while (n <= length) {
      const placed = !!st.covered[n], col = placed ? colourFn(n) : null;
      let m = n;
      while (m + 1 <= length && !!st.covered[m + 1] === placed && (!placed || sameColour(colourFn(m + 1), col))) m++;
      if (placed) kids.push(rect(n - 1, m - n + 1, 6, 16, css(col)));
      else kids.push(rect(n - 1, m - n + 1, 13, 2, '#8a94a8', { opacity: '0.7' }));
      n = m + 1;
    }
    if (st.res508) kids.push(rect(502, 12, 2, 24, css(ROSE), { 'data-mark': '508' }));
    else kids.push(rect(502, 12, 3, 22, 'none', { stroke: css(ROSE), 'stroke-width': '2', 'data-mark': '508-absent' }));
    el.ruler.replaceChildren(...kids);
  }

  // ---- drawing ----
  function size() {
    dpr = Math.min(window.devicePixelRatio || 1, 2);
    const r = canvas.getBoundingClientRect();
    seenW = r.width; seenH = r.height;                                       // the size this frame was drawn at; the resize observer compares against it
    w = Math.max(1, Math.round(r.width * dpr)); h = Math.max(1, Math.round(r.height * dpr));
    if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
  }
  function mvp() {
    const dist = (radius / Math.tan(FOV / 2)) * 1.15 / S.zoom;
    const sp = Math.sin(S.phi), eye = [dist * sp * Math.cos(S.theta), dist * sp * Math.sin(S.theta), dist * Math.cos(S.phi)];
    const view = lookAt(eye, [0, 0, 1]);
    return { m: multiply(perspective(FOV, w / h, dist * 0.1, dist + radius * 3), view), view, dist };
  }
  function walkStrength() {              // the walk's emphasis holds for a moment, then eases away, so no structure stays darker than another for no reason in the data
    if (S.walkAt === 0) return 0;
    const idle = performance.now() - S.walkAt;
    return Math.max(0, Math.min(1, 1 - (idle - 1800) / 900));
  }
  function drawStructure(idx, m, view, dist, o) {
    // One structure, complete, with its own depth: near parts of it hide far parts of it. The caller clears depth between structures, so a faded or ghosted
    // structure never hides the one drawn after it. A structure at partial alpha is drawn depth-first and colour-second, so each pixel shows only the nearest
    // piece of it once, rather than stacking translucent layers.
    const g = gpu.structures[idx], u = gpu.tubeU;
    gl.useProgram(gpu.tube);
    gl.uniformMatrix4fv(u.uMvp, false, m); gl.uniformMatrix4fv(u.uView, false, view); gl.uniform1f(u.uDist, dist); gl.uniform1f(u.uRad, radius);
    gl.uniform1f(u.uProgress, o.progress); gl.uniform1f(u.uWalk, o.walk); gl.uniform1f(u.uWalkK, walkStrength()); gl.uniform1f(u.uComet, 30);
    gl.uniform1f(u.uGrey, o.grey ? 1 : 0); gl.uniform3fv(u.uGreyColour, [0.5, 0.52, 0.56]); gl.uniform1f(u.uAlpha, o.alpha);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, gpu.palette); gl.uniform1i(u.uPal, 0); gl.uniform1i(u.uPalMax, gpu.paletteN - 1);
    gl.bindVertexArray(g.vao);
    gl.enable(gl.DEPTH_TEST); gl.depthMask(true);
    draws.depthOn = gl.isEnabled(gl.DEPTH_TEST); draws.tube++;
    if (o.alpha < 1) {
      gl.colorMask(false, false, false, false); gl.depthFunc(gl.LESS);
      gl.drawElements(gl.TRIANGLES, g.count, g.type, 0);
      gl.colorMask(true, true, true, true); gl.depthFunc(gl.LEQUAL);
      gl.drawElements(gl.TRIANGLES, g.count, g.type, 0);
      gl.depthFunc(gl.LESS);
    } else {
      gl.drawElements(gl.TRIANGLES, g.count, g.type, 0);
    }
    gl.bindVertexArray(null);
    if (g.dotCount) {
      const p = gpu.pointU;
      gl.useProgram(gpu.point);
      gl.uniformMatrix4fv(p.uMvp, false, m); gl.uniform1f(p.uUseAttr, 1); gl.uniform1f(p.uSize, 3 * dpr); gl.uniform3fv(p.uColour, GREY); gl.uniform1f(p.uMode, 0); gl.uniform1f(p.uAlpha, 0.55 * o.alpha);
      gl.depthMask(false);                                                   // the dots sit in the scene (the tube can hide them) but hide nothing themselves
      gl.bindVertexArray(g.dots);
      draws.dots++;
      gl.drawArrays(gl.POINTS, 0, g.dotCount);
      gl.bindVertexArray(null);
      gl.depthMask(true);
    }
  }
  function marker(m, pos, css, colour, mode, alpha = 1) {
    const u = gpu.pointU, range = gl.getParameter(gl.ALIASED_POINT_SIZE_RANGE);
    gl.useProgram(gpu.point);
    gl.uniformMatrix4fv(u.uMvp, false, m); gl.uniform1f(u.uUseAttr, 0); gl.uniform3fv(u.uPos, pos); gl.uniform1f(u.uSize, Math.min(range[1], css * dpr));
    gl.uniform3fv(u.uColour, colour); gl.uniform1f(u.uMode, mode); gl.uniform1f(u.uAlpha, alpha);
    gl.drawArrays(gl.POINTS, 0, 1);
  }
  function walkResidue() {
    if (S.intro) return Math.max(1, Math.floor(S.intro.progress));
    return S.walk;
  }
  function draw() {
    if (!gpu || lost) return;
    size();
    gl.viewport(0, 0, w, h);
    gl.clearColor(0, 0, 0, 0); gl.clearDepth(1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    const { m, view, dist } = mvp();
    const st = data.structures[S.i];
    const progress = S.intro ? S.intro.progress : 1e9;
    const walk = S.walk !== null && !S.intro ? S.walk : -1;
    if (S.overlay && S.i !== 0) { drawStructure(0, m, view, dist, { alpha: 0.35, grey: true, progress: 1e9, walk: -1 }); gl.clear(gl.DEPTH_BUFFER_BIT); }
    if (S.fade) {
      drawStructure(S.fade.from, m, view, dist, { alpha: 1 - S.fade.t, progress, walk: -1 });
      gl.clear(gl.DEPTH_BUFFER_BIT);
      drawStructure(S.i, m, view, dist, { alpha: S.fade.t, progress, walk });
    } else {
      drawStructure(S.i, m, view, dist, { alpha: 1, progress, walk });
    }
    gl.disable(gl.DEPTH_TEST);                                               // markers are annotations: always visible, drawn over the shape
    const ready = !S.intro || S.intro.progress >= 508;
    const fadeA = S.fade ? S.fade.t : 1;
    if (ready) {
      if (st.res508) { marker(m, st.res508.xyz, 34, ROSE, 1, fadeA); marker(m, st.res508.xyz, 16, ROSE, 0, fadeA); }
      else if (st.at.has(507) && st.at.has(509)) {                       // where 508 would be: midway between its neighbours, drawn as an empty ring
        const a = st.at.get(507), b = st.at.get(509);
        draws.absentRingMode = 3;                                             // 3 is the dashed ring: nothing is placed here
        marker(m, [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2], 30, ROSE, 3, fadeA);
        if (S.overlay && data.structures[0].res508) marker(m, data.structures[0].res508.xyz, 22, GREY, 1, 0.8);     // where 6MSM places it, as a grey ring
      }
    }
    for (const fn of layerFns.values()) {
      for (const mk of fn(st, S.i) || []) marker(m, mk.pos, mk.size, mk.colour, mk.mode ?? 0, mk.alpha ?? 1);
    }
    const wr = walkResidue();
    if (wr !== null && st.at.has(wr) && !(st.res508 && wr === 508 && !S.intro)) { marker(m, st.at.get(wr), 38, WHITE, 2); marker(m, st.at.get(wr), 11, WHITE, 0); }
    else if (wr === 508 && st.res508 && !S.intro) marker(m, st.res508.xyz, 9, WHITE, 0);
    canvas.dataset.frames = String(++frames);
    canvas.dataset.draws = JSON.stringify(draws);
    lastDrawT = performance.now();
    dirty = false;
  }
  function frame(t) {
    raf = 0;
    if (document.hidden || offscreen) return;
    const dt = lastT ? Math.min(0.05, (t - lastT) / 1000) : 0;
    lastT = t;
    let animating = false;
    if (S.intro) {
      S.intro.progress = Math.min(length + 1, ((t - S.intro.start) / INTRO_MS) * length);
      if (S.intro.progress >= length + 1) { S.intro = null; publish(); }
      animating = true; dirty = true;
    }
    if (S.walk !== null && S.walkAt && walkStrength() > 0 && walkStrength() < 1) { animating = true; dirty = true; }
    else if (S.walkAt && walkStrength() === 0 && S.walk !== null && !S._walkFaded) { S._walkFaded = true; dirty = true; }
    if (S.fade) {
      S.fade.t = Math.min(1, (t - S.fade.start) / FADE_MS);
      if (S.fade.t >= 1) { S.fade = null; if (el.note) el.note.textContent = ''; publish(); }
      animating = true; dirty = true;
    }
    if (S.playing && !S.intro) { S.theta += dt * TURN_RATE; syncSliders(); }
    // while only the slow turn is running, draw at most about 30 times a second
    if (dirty || animating || !S.playing || t - lastDrawT >= 32) draw();
    if (S.playing || animating) schedule(); else lastT = 0;
  }
  function schedule() { if (!raf && !document.hidden && !offscreen && !lost) raf = requestAnimationFrame(frame); }
  function request() { dirty = true; schedule(); }

  // ---- controls ----
  function syncSliders() {
    if (el.turn) el.turn.value = String(Math.round(((S.theta % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI) * 180 / Math.PI));
    if (el.tilt) el.tilt.value = String(Math.round(S.phi * 100) / 100);
  }
  function syncPlay() { if (el.play) el.play.setAttribute('aria-pressed', S.playing ? 'true' : 'false'); }
  function select(i, animate = false) {
    const from = S.i;
    S.i = i;
    if (animate && i !== from && !reduceQuery.matches) {
      S.fade = { from, t: 0, start: performance.now() };
      if (el.note) el.note.textContent = 'Cross-fade between two separate structures, ' + data.structures[from].id + ' and ' + data.structures[i].id + ', not a movement.';
    }
    if (el.picks) for (const b of el.picks.querySelectorAll('button')) b.setAttribute('aria-pressed', String(Number(b.dataset.index) === i));
    if (el.overlay) { el.overlay.disabled = i === 0; el.overlay.setAttribute('aria-pressed', String(S.overlay && i !== 0)); }
    updateReadout(); request();
  }
  function setOverlay(on_) { S.overlay = on_; if (el.overlay) el.overlay.setAttribute('aria-pressed', String(on_ && S.i !== 0)); updateReadout(); request(); }
  function setWalk(n) { S.walk = n; S._walkFaded = false; S.walkAt = n === null ? 0 : performance.now(); if (n !== null && el.walk) el.walk.value = String(n); updateReadout(); request(); if (n !== null) { setTimeout(request, 1850); setTimeout(request, 2800); } }
  function view(dTheta, dPhi) { S.theta += dTheta; S.phi = Math.max(0.2, Math.min(Math.PI - 0.2, S.phi + dPhi)); S.viewTouched = true; syncSliders(); request(); publish(); }
  function applyHash() {
    if (HASH_MODE !== 'url') return;
    const q = Object.fromEntries(location.hash.replace(/^#/, '').split('&').filter(Boolean).map(kv => kv.split('=')));
    const num = (k, lo, hi) => { const v = Number(q[k]); return q[k] !== undefined && Number.isFinite(v) && v >= lo && v <= hi ? v : null; };
    if (q.s !== undefined && byId[q.s] !== undefined) S.i = byId[q.s];
    const r = num('r', 1, length); if (r !== null && Number.isInteger(r)) { S.walk = r; if (el.walk) el.walk.value = String(r); }
    if (q.o === '1' && S.i !== 0) S.overlay = true;
    const y = num('y', 0, 360), p = num('p', 11, 169), z = num('z', 0.5, 4);
    if (y !== null && p !== null && z !== null) { S.theta = y * Math.PI / 180; S.phi = p * Math.PI / 180; S.zoom = z; S.viewTouched = true; S.playing = false; }
    if (q.l) for (const name of q.l.split(',')) if (/^[a-z]{1,16}$/.test(name)) S.layers.add(name);
  }
  function reset() { Object.assign(S, home); S.viewTouched = false; syncSliders(); request(); publish(); }
  function zoomBy(f) { S.zoom = Math.max(0.5, Math.min(4, S.zoom * f)); S.viewTouched = true; request(); publish(); }

  if (el.picks) {
    data.structures.forEach((st, i) => {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'cftr-btn'; b.dataset.index = String(i); b.textContent = st.id; b.setAttribute('aria-pressed', i === 0 ? 'true' : 'false');
      b.addEventListener('click', () => select(i, true));
      el.picks.append(b);
    });
  }
  function legendRow(colour, text, kind) {
    const li = document.createElement('li'), sw = document.createElement('span');
    sw.className = 'cftr-swatch' + (kind ? ' is-' + kind : '');
    if (kind === 'dotted') { sw.style.background = 'transparent'; sw.style.borderTop = '2px dotted ' + colour; } else sw.style.background = colour;
    sw.setAttribute('aria-hidden', 'true');
    li.append(sw, text);
    return li;
  }
  if (el.legend) {
    el.legend.replaceChildren(...parts.map(p => {
      const li = document.createElement('li'), sw = document.createElement('span');
      sw.className = 'cftr-swatch'; sw.style.background = 'rgb(' + p.colour.map(v => Math.round(v * 255)).join(',') + ')'; sw.setAttribute('aria-hidden', 'true');
      li.append(sw, p.label + ', residues ' + p.start + ' to ' + p.end);
      return li;
    }).concat([legendRow('#8c99a8', 'Grey: placed, outside the four marked parts'), legendRow('#8c99a8', 'Dotted grey: not placed in this structure', 'dotted'),
               legendRow('#c9d1dc', 'Flat band: helix or strand, as UniProt annotates', 'band')]));
  }
  if (el.walk) { el.walk.min = '1'; el.walk.max = String(length); el.walk.value = '508'; }
  on(el.walk, 'input', () => setWalk(Number(el.walk.value)));
  on($('cftr-walk-clear'), 'click', () => setWalk(null));
  on($('cftr-walk-508'), 'click', () => setWalk(508));
  on(el.overlay, 'click', () => setOverlay(!S.overlay));
  on(el.play, 'click', () => { S.playing = !S.playing; syncPlay(); publish(); request(); });
  on(el.turn, 'input', () => { S.theta = Number(el.turn.value) * Math.PI / 180; S.viewTouched = true; request(); publish(); });
  on(el.tilt, 'input', () => { S.phi = Number(el.tilt.value); S.viewTouched = true; request(); publish(); });
  on($('cftr-zoom-in'), 'click', () => zoomBy(1.25));
  on($('cftr-zoom-out'), 'click', () => zoomBy(1 / 1.25));
  on($('cftr-reset'), 'click', reset);
  reduceQuery.addEventListener('change', () => { S.playing = !reduceQuery.matches && !turnOff; if (reduceQuery.matches) { S.intro = null; S.fade = null; } syncPlay(); publish(); request(); });

  let drag = null;
  canvas.addEventListener('pointerdown', (e) => { drag = { x: e.clientX, y: e.clientY, touch: e.pointerType === 'touch' }; try { canvas.setPointerCapture(e.pointerId); } catch (x) { /* not capturable */ } });
  canvas.addEventListener('pointermove', (e) => {
    if (!drag) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag.x = e.clientX; drag.y = e.clientY;
    view(dx * 0.01, drag.touch ? 0 : dy * 0.01);                          // a touch drag turns; tilting is on the slider so a vertical swipe can still scroll the page
  });
  const endDrag = (e) => { drag = null; try { canvas.releasePointerCapture(e.pointerId); } catch (x) { /* already released */ } };
  canvas.addEventListener('pointerup', endDrag); canvas.addEventListener('pointercancel', endDrag);
  canvas.addEventListener('keydown', (e) => {
    const k = e.key, step = 0.12;
    if (k === 'ArrowLeft') view(-step, 0); else if (k === 'ArrowRight') view(step, 0); else if (k === 'ArrowUp') view(0, -step); else if (k === 'ArrowDown') view(0, step);
    else if (k === '+' || k === '=') zoomBy(1.25); else if (k === '-') zoomBy(1 / 1.25);
    else if (k === 'r' || k === 'R') reset();
    else if (/^[1-9]$/.test(k) && Number(k) <= data.structures.length) select(Number(k) - 1, true);
    else return;
    e.preventDefault();
  });
  canvas.addEventListener('wheel', (e) => { if (document.activeElement !== canvas) return; e.preventDefault(); zoomBy(e.deltaY < 0 ? 1.1 : 1 / 1.1); }, { passive: false });

  for (const li of document.querySelectorAll('[data-step]')) {
    const btn = li.querySelector('button');
    if (!btn) continue;
    btn.hidden = false;
    btn.addEventListener('click', () => {
      for (const o of document.querySelectorAll('[data-step]')) o.removeAttribute('aria-current');
      li.setAttribute('aria-current', 'step');
      S.overlay = li.dataset.overlay === '1';
      S.walk = li.dataset.walk ? Number(li.dataset.walk) : null;
      if (S.walk !== null && el.walk) el.walk.value = String(S.walk);
      select(byId[li.dataset.structure] ?? 0, true);
    });
  }

  const publishMesh = () => { canvas.dataset.mesh = JSON.stringify({ triangles: gpu.triangles, dots: gpu.structures.map(g => g.dotCount), buildMs: Math.round(gpu.buildMs * 10) / 10, meshMs: Math.round(gpu.meshMs * 10) / 10 }); };
  canvas.addEventListener('webglcontextlost', (e) => { e.preventDefault(); lost = true; if (raf) cancelAnimationFrame(raf); raf = 0; });
  canvas.addEventListener('webglcontextrestored', () => { gpu = build(gl, data, colourFn, sec); publishMesh(); lost = false; request(); });
  document.addEventListener('visibilitychange', () => { if (document.hidden) { if (raf) cancelAnimationFrame(raf); raf = 0; lastT = 0; } else request(); });
  if (window.IntersectionObserver) new IntersectionObserver((entries) => {     // off the screen: stop the loop; back on it: one redraw, and the turn resumes if it was on
    const seen = entries[entries.length - 1].isIntersecting;
    if (seen === !offscreen) return;
    offscreen = !seen;
    if (offscreen) { if (raf) cancelAnimationFrame(raf); raf = 0; lastT = 0; } else request();
    publish();
  }, { threshold: 0 }).observe(canvas);
  // a size change asks for a redraw; the observer's first report (and any report of the same size) does not, so a still picture stays still
  const resized = () => { const r = canvas.getBoundingClientRect(); return r.width !== seenW || r.height !== seenH; };
  if (window.ResizeObserver) new ResizeObserver(() => { if (resized()) request(); }).observe(canvas); else window.addEventListener('resize', () => { if (resized()) request(); });

  try { gpu = build(gl, data, colourFn, sec); } catch (e) { return fail('The 3D view could not start (' + e.message + ').'); }
  publishMesh();
  if (root) root.hidden = false;
  if (fallback) fallback.hidden = true;
  for (const n of document.querySelectorAll('[data-needs-gl]')) n.hidden = false;
  if (hero && !reduceQuery.matches) S.intro = { start: performance.now(), progress: 0 };      // the chain draws itself once, in residue order
  applyHash();
  syncSliders(); syncPlay(); select(S.i);
  draw();                                                                  // one frame straight away, so a visitor who asked for less motion still sees the shape
  schedule();
  return {
    data, parts, canvas, length,
    state: snapshot,
    onChange: (fn) => { listeners.push(fn); },
    select: (id) => { if (byId[id] !== undefined) select(byId[id], true); },
    setWalk, setOverlay,
    setView: (v) => { if (v.theta !== undefined) S.theta = v.theta; if (v.phi !== undefined) S.phi = v.phi; if (v.zoom !== undefined) S.zoom = v.zoom; S.viewTouched = true; syncSliders(); request(); publish(); },
    setColour: (fn) => { colourFn = fn || ((n) => colourOf(parts, n)); if (gpu) recolour(gl, gpu, colourFn); updateReadout(); request(); },
    setLayer: (name, on_, fn) => { if (on_) { S.layers.add(name); if (fn) layerFns.set(name, fn); } else { S.layers.delete(name); layerFns.delete(name); } request(); publish(); },
    layerNames: () => [...S.layers],
    request,
    redrawNow: () => { dirty = true; draw(); },
    structureIndex: (id) => byId[id],
    encodeHash,
  };
}

// The pure geometry, exported so tests can check it with synthetic input (the tube passes through every measured point, a helix is wider than a coil, the counts add up).
export const internals = { tubeMesh, splinePositions, catmullRom, secondaryByResidue, crowding, SHAPES, SUBDIV, SIDES, SEC_CODE };
export const ready = start();
