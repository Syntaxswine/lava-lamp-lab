// render.js — one WebGL2 pass. Everything is analytic or volumetric: the glass
// globe as a surface of revolution, the metal base and cap as signed-distance
// cones, and the wax as a raymarched 3-D field.
//
// The wax is lit the way a real lamp lights it: one warm source directly under
// the globe, shining up through the fluid. Because the light travels along +y,
// the optical depth reaching any voxel is a prefix sum down its own column, and
// volume.js has already baked that into the alpha channel. So the blobs shade
// themselves and shadow the blobs above them for the cost of the fetch we were
// making anyway.

import { GEO, WAX } from './params.js';
import { Volume, TMIN, TMAX } from './volume.js';

const VS = `#version 300 es
void main() {
  vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2) * 2.0 - 1.0;
  gl_Position = vec4(p, 0.0, 1.0);
}`;

const FS = `#version 300 es
precision highp float;
precision highp sampler3D;

uniform vec2  uRes;
uniform vec3  uCam;
uniform mat3  uBasis;          // right, up, forward
uniform float uTan;
uniform sampler3D uVol;
uniform vec3  uBoxMin, uBoxMax;
uniform float uProf[64];
uniform float uH, uRmax;
uniform vec3  uWax, uWaxDeep, uFluid, uBulb;
uniform float uSteps, uExposure, uGlow;
uniform int   uView;           // 0 lamp, 1 temperature, 2 rise/sink
uniform float uTcross;         // crossover, normalised into TMIN..TMAX
uniform float uTfluidTop, uTfluidBase;
uniform float uTmin, uTmax;
uniform float uBulbOn;
uniform float uTmelt, uDTmelt;   // to derive solid fraction from temperature
out vec4 frag;

const float PI = 3.14159265;

float profileR(float y) {
  float t = clamp(y / uH, 0.0, 1.0) * 63.0;
  int i = int(floor(t));
  int j = min(i + 1, 63);
  return mix(uProf[i], uProf[j], t - float(i));
}

// distance from the axis, minus the glass radius at that height
float insideGlobe(vec3 p) {
  if (p.y < 0.0 || p.y > uH) return -1.0;
  return profileR(p.y) - length(p.xz);
}

// outward normal of the surface of revolution
vec3 globeNormal(vec3 p) {
  float e = 0.002;
  float dr = (profileR(p.y + e) - profileR(p.y - e)) / (2.0 * e);
  vec2 radial = normalize(p.xz + vec2(1e-6));
  return normalize(vec3(radial.x, -dr, radial.y));
}

// ---- ray/bounding-cylinder --------------------------------------------------
bool cylinder(vec3 ro, vec3 rd, float R, float y0, float y1, out float t0, out float t1) {
  float a = dot(rd.xz, rd.xz);
  float b = dot(ro.xz, rd.xz);
  float c = dot(ro.xz, ro.xz) - R * R;
  float h = b * b - a * c;
  if (h < 0.0) return false;
  h = sqrt(h);
  t0 = (-b - h) / a;
  t1 = (-b + h) / a;
  float ta = (y0 - ro.y) / rd.y;
  float tb = (y1 - ro.y) / rd.y;
  if (ta > tb) { float s = ta; ta = tb; tb = s; }
  t0 = max(t0, ta); t1 = min(t1, tb);
  return t1 > max(t0, 0.0);
}

// ---- metal fittings, as signed distance ------------------------------------
float sdCappedCone(vec3 p, float h, float r1, float r2) {
  vec2 q = vec2(length(p.xz), p.y);
  vec2 k1 = vec2(r2, h);
  vec2 k2 = vec2(r2 - r1, 2.0 * h);
  vec2 ca = vec2(q.x - min(q.x, (q.y < 0.0) ? r1 : r2), abs(q.y) - h);
  vec2 cb = q - k1 + k2 * clamp(dot(k1 - q, k2) / dot(k2, k2), 0.0, 1.0);
  float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
  return s * sqrt(min(dot(ca, ca), dot(cb, cb)));
}

float ribs(vec3 p) {
  float a = atan(p.z, p.x);
  return 0.0007 * sin(a * 22.0) * smoothstep(0.0, 0.015, abs(p.y));
}

float metal(vec3 p) {
  // proportions follow the globe: a base cone flaring from the foot, and a
  // matching cap flaring up from the neck
  vec3 b = p - vec3(0.0, -0.040, 0.0);
  float base = sdCappedCone(b, 0.040, 0.046, 0.0280) - ribs(p);
  vec3 c = p - vec3(0.0, uH + 0.022, 0.0);
  float cap = sdCappedCone(c, 0.022, 0.0128, 0.0215) - ribs(p);
  return min(base, cap);
}

vec3 metalNormal(vec3 p) {
  vec2 e = vec2(0.0004, 0.0);
  return normalize(vec3(
    metal(p + e.xyy) - metal(p - e.xyy),
    metal(p + e.yxy) - metal(p - e.yxy),
    metal(p + e.yyx) - metal(p - e.yyx)));
}

vec3 metalShade(vec3 p, vec3 rd) {
  vec3 n = metalNormal(p);
  vec3 l = normalize(vec3(-0.5, 0.85, 0.4));
  float diff = max(0.0, dot(n, l));
  float rim = pow(1.0 - abs(dot(n, -rd)), 3.0);
  vec3 base = vec3(0.335, 0.245, 0.148);           // aged brass, not bone
  vec3 mcol = base * (0.16 + 0.72 * diff) + vec3(1.0, 0.88, 0.68) * pow(diff, 44.0) * 0.20;
  mcol += vec3(0.55, 0.60, 0.72) * rim * 0.10;
  return mcol;
}

// ---- palettes ---------------------------------------------------------------
vec3 heat(float t) {                       // temperature view
  t = clamp(t, 0.0, 1.0);
  return vec3(smoothstep(0.35, 1.0, t), smoothstep(0.05, 0.75, t) * 0.85,
              1.0 - smoothstep(0.15, 0.7, t));
}

void main() {
  vec2 uv = (gl_FragCoord.xy / uRes) * 2.0 - 1.0;
  uv.x *= uRes.x / uRes.y;
  vec3 rd = normalize(uBasis * vec3(uv * uTan, 1.0));
  vec3 ro = uCam;

  // background: a dark room with the lamp glowing into it
  float hor = smoothstep(-0.35, 0.5, rd.y);
  vec3 col = mix(vec3(0.0075, 0.0100, 0.0150), vec3(0.0018, 0.0026, 0.0044), hor);
  vec3 toBulb = vec3(0.0, -0.02, 0.0) - ro;
  float bulbAlign = max(0.0, dot(rd, normalize(toBulb)));
  col += uBulb * uBulbOn * 0.10 * pow(bulbAlign, 22.0);

  // ---- metal fittings ------------------------------------------------------
  float tMetal = 1e9;
  {
    float t = 0.05;
    for (int i = 0; i < 72; i++) {
      vec3 p = ro + rd * t;
      float d = metal(p);
      if (d < 0.0002) { tMetal = t; break; }
      t += max(d, 0.0006);
      if (t > 3.0) break;
    }
  }

  // ---- glass + volume -----------------------------------------------------
  // Depth matters here. The base plate and the coil sit INSIDE the globe, so a
  // ray can enter the glass, cross wax, and only then hit metal. Drawing the
  // metal after the volume and letting it win painted the plate over the pool in
  // front of it -- a bright bone-white disc through several centimetres of wax.
  // The metal is the BACKDROP the volume composites over, unless it is nearer
  // than the glass.
  float t0, t1;
  bool hitGlobe = cylinder(ro, rd, uRmax + 0.004, 0.0, uH, t0, t1);
  bool metalFront = tMetal < 1e8 && (!hitGlobe || tMetal <= max(t0, 0.0));
  vec3 backdrop = col;
  // Shade the metal whenever it is hit, not only when it is behind the glass.
  // Guarding this on !metalFront left every ray that met the base or the cap
  // directly -- which is every ray that sees them at all -- compositing the
  // BACKGROUND colour, and the fittings simply vanished from the picture.
  if (tMetal < 1e8) {
    vec3 pm = ro + rd * tMetal;
    backdrop = metalShade(pm, rd);
    backdrop += uBulb * uBulbOn * 0.07 * exp(-34.0 * abs(pm.y + 0.005));
  }
  vec3 lamp = vec3(0.0);
  float residual = 1.0;
  vec3 glassSpec = vec3(0.0);

  if (hitGlobe && !metalFront) {
    t0 = max(t0, 0.0);
    if (tMetal < 1e8) t1 = min(t1, tMetal);
    float steps = uSteps;
    float ds = (t1 - t0) / steps;
    vec3 inv = 1.0 / (uBoxMax - uBoxMin);
    vec3 acc = vec3(0.0);
    float trans = 1.0;
    bool haveSurface = false;
    vec3 sNormal = vec3(0.0);
    vec3 sPos = vec3(0.0);
    float entered = -1.0;
    float fluidPath = 0.0;
    float tOutside = t0;

    for (float i = 0.5; i < steps; i += 1.0) {
      float t = t0 + i * ds;
      vec3 p = ro + rd * t;
      if (insideGlobe(p) < 0.0) { tOutside = t; continue; }
      if (entered < 0.0) {
        // Find the glass surface, do not settle for the first sample that
        // happened to land inside it. At 96 steps across the bore that sample can
        // be a centimetre deep, and the "surface normal" taken there points
        // somewhere between the glass and the axis -- which put a Fresnel sheen
        // across the whole bottle and made the fluid look like milk.
        float lo = tOutside, hi = t;
        for (int b = 0; b < 6; b++) {
          float mid = 0.5 * (lo + hi);
          if (insideGlobe(ro + rd * mid) < 0.0) lo = mid; else hi = mid;
        }
        entered = hi;
      }
      fluidPath += ds;

      vec4 s = texture(uVol, (p - uBoxMin) * inv);
      float wax = s.r;
      if (wax > 0.004) {
        // gradient normal at the first substantial crossing, for the highlight
        if (!haveSurface && wax > 0.16) {
          vec3 e = (uBoxMax - uBoxMin) / vec3(textureSize(uVol, 0));
          float dx = texture(uVol, (p + vec3(e.x, 0, 0) - uBoxMin) * inv).r
                   - texture(uVol, (p - vec3(e.x, 0, 0) - uBoxMin) * inv).r;
          float dy = texture(uVol, (p + vec3(0, e.y, 0) - uBoxMin) * inv).r
                   - texture(uVol, (p - vec3(0, e.y, 0) - uBoxMin) * inv).r;
          float dz = texture(uVol, (p + vec3(0, 0, e.z) - uBoxMin) * inv).r
                   - texture(uVol, (p - vec3(0, 0, e.z) - uBoxMin) * inv).r;
          vec3 g = vec3(dx, dy, dz);
          if (dot(g, g) > 1e-9) { sNormal = -normalize(g); sPos = p; haveSurface = true; }
        }

        float shade = s.a;                        // light that reached this voxel
        float Tn = s.g;
        // solid fraction, from the same clamp wax.js uses -- not stored
        float Tc = uTmin + Tn * (uTmax - uTmin);
        float solid = clamp((uTmelt + uDTmelt - Tc) / (2.0 * uDTmelt), 0.0, 1.0);

        vec3 body;
        if (uView == 1) {
          body = heat(Tn);
        } else if (uView == 2) {
          // does this parcel rise or sink where it is? compare its temperature
          // with the fluid it displaces at its own height
          float fl = mix(uTfluidBase, uTfluidTop, clamp(p.y / uH, 0.0, 1.0));
          float lift = (Tn - uTcross) - 0.45 * (fl - uTcross);
          body = mix(vec3(0.15, 0.42, 0.95), vec3(1.0, 0.45, 0.12),
                     smoothstep(-0.06, 0.06, lift));
        } else {
          // dyed wax: deep colour in the shadowed core, bright where the bulb
          // shines through, and pale and matte where it has not melted
          body = mix(uWaxDeep, uWax, shade);
          body += uBulb * uBulb * shade * shade * 0.55;
          body = mix(body, mix(uWax, vec3(0.93, 0.90, 0.86), 0.72), solid * 0.85);
        }

        // Extinction of dyed wax. 26 per metre left a 2 cm blob only 22%
        // opaque and the lamp looked like smoke; 110 per metre puts it at 70%,
        // which is what a real blob looks like held against a bulb.
        float sigma = 110.0 * wax;
        float a = 1.0 - exp(-sigma * ds);
        // A blob deep in another blob's shadow is not black: light diffuses
        // through wax rather than stopping at it, which is why the whole globe
        // of a running lamp glows and not just the bottom blob.
        acc += trans * a * body * (0.30 + 0.80 * (uView == 0 ? shade : 1.0));
        if (uView == 0) acc += trans * a * uWaxDeep * 0.45;
        trans *= 1.0 - a;
        if (trans < 0.008) break;
      }
    }

    if (entered >= 0.0) {
      // the fluid itself: faint tint, and it carries the bulb glow
      vec3 pEnter = ro + rd * entered;
      float ft = 1.0 - exp(-0.55 * fluidPath);
      float lampGlow = uBulbOn * exp(-6.5 * max(0.0, pEnter.y));
      vec3 fluidCol = uFluid * (0.030 + 0.55 * lampGlow) + uBulb * lampGlow * 0.30;
      acc += trans * ft * fluidCol;
      trans *= 1.0 - ft;
      residual = trans;

      // glass: Fresnel rim on the entry surface
      vec3 n = globeNormal(pEnter);
      float fres = pow(1.0 - abs(dot(n, -rd)), 4.0);
      glassSpec = vec3(0.52, 0.64, 0.80) * fres * 0.34;
      vec3 hv = normalize(normalize(vec3(-0.45, 0.75, 0.5)) - rd);
      glassSpec += vec3(1.0) * pow(max(0.0, dot(n, hv)), 110.0) * 0.35;

      // specular on the wax surface, from the bulb below and the room light
      if (haveSurface && uView == 0) {
        vec3 ldir = normalize(vec3(0.0, -0.10, 0.0) - sPos);
        vec3 hb = normalize(ldir - rd);
        float sp = pow(max(0.0, dot(sNormal, hb)), 42.0);
        acc += uBulb * sp * 0.55 * uBulbOn;
        vec3 hr = normalize(normalize(vec3(-0.5, 0.8, 0.35)) - rd);
        acc += vec3(0.9, 0.95, 1.0) * pow(max(0.0, dot(sNormal, hr)), 60.0) * 0.35;
      }
    }
    lamp = acc;
    residual = trans;
  }

  // composite: whatever the volume did not absorb shows the backdrop through it
  if (metalFront) {
    col = backdrop;
  } else if (hitGlobe) {
    col = backdrop * residual + lamp + glassSpec;
  } else {
    col = backdrop;
  }

  // glow around the lamp, cheap: the emissive we already accumulated
  col += lamp * lamp * uGlow;

  col *= uExposure;
  col = col / (1.0 + col);                       // tone map
  col = pow(col, vec3(1.0 / 2.2));
  frag = vec4(col, 1.0);
}`;

export const PALETTES = {
  classic:   { wax: [1.00, 0.20, 0.10], deep: [0.32, 0.03, 0.03], fluid: [0.95, 0.72, 0.20], bulb: [1.00, 0.62, 0.26] },
  tangerine: { wax: [1.00, 0.44, 0.06], deep: [0.30, 0.07, 0.01], fluid: [0.10, 0.42, 0.55], bulb: [1.00, 0.70, 0.34] },
  violet:    { wax: [1.00, 0.84, 0.16], deep: [0.34, 0.20, 0.02], fluid: [0.36, 0.12, 0.62], bulb: [1.00, 0.80, 0.42] },
  glacier:   { wax: [0.94, 0.97, 1.00], deep: [0.24, 0.34, 0.46], fluid: [0.06, 0.35, 0.68], bulb: [0.86, 0.94, 1.00] },
  absinthe:  { wax: [0.62, 1.00, 0.28], deep: [0.10, 0.26, 0.05], fluid: [0.62, 0.14, 0.34], bulb: [0.92, 1.00, 0.62] },
};

export class Renderer {
  constructor(canvas) {
    this.canvas = canvas;
    const gl = canvas.getContext('webgl2', {
      antialias: false, alpha: false, powerPreference: 'high-performance',
    });
    if (!gl) throw new Error('WebGL2 unavailable');
    this.gl = gl;
    this.volume = new Volume();
    this.scale = 0.72;
    this.steps = 96;
    this.palette = PALETTES.classic;
    this.view = 0;
    this.exposure = 1.25;
    this.glow = 0.16;
    this.bulbOn = 1;
    this.orbit = { az: 0.6, el: 0.12, dist: 0.60 };
    this.buildProgram();
    this.buildTexture();
  }

  buildProgram() {
    const gl = this.gl;
    const sh = (type, src) => {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
        throw new Error(gl.getShaderInfoLog(s) || 'shader compile failed');
      }
      return s;
    };
    const p = gl.createProgram();
    gl.attachShader(p, sh(gl.VERTEX_SHADER, VS));
    gl.attachShader(p, sh(gl.FRAGMENT_SHADER, FS));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
      throw new Error(gl.getProgramInfoLog(p) || 'link failed');
    }
    this.prog = p;
    gl.useProgram(p);
    this.u = {};
    const names = ['uRes', 'uCam', 'uBasis', 'uTan', 'uVol', 'uBoxMin', 'uBoxMax',
      'uProf', 'uH', 'uRmax', 'uWax', 'uWaxDeep', 'uFluid', 'uBulb', 'uSteps',
      'uExposure', 'uGlow', 'uView', 'uTcross', 'uTfluidTop', 'uTfluidBase', 'uBulbOn',
      'uTmelt', 'uDTmelt', 'uTmin', 'uTmax'];
    for (const n of names) this.u[n] = gl.getUniformLocation(p, n);
    gl.uniform1fv(this.u.uProf, this.volume.profile);
    gl.uniform1f(this.u.uH, GEO.H);
    gl.uniform1f(this.u.uRmax, this.volume.rmax);
    gl.uniform3fv(this.u.uBoxMin, this.volume.min);
    gl.uniform3fv(this.u.uBoxMax, this.volume.max);
    gl.uniform1i(this.u.uVol, 0);
    gl.uniform1f(this.u.uTmelt, WAX.Tmelt);
    gl.uniform1f(this.u.uDTmelt, WAX.dTmelt);
    gl.uniform1f(this.u.uTmin, TMIN);
    gl.uniform1f(this.u.uTmax, TMAX);
    this.vao = gl.createVertexArray();
  }

  buildTexture() {
    const gl = this.gl, v = this.volume;
    this.tex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_3D, this.tex);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_R, gl.CLAMP_TO_EDGE);
    gl.texImage3D(gl.TEXTURE_3D, 0, gl.RGBA8, v.nx, v.ny, v.nz, 0,
      gl.RGBA, gl.UNSIGNED_BYTE, null);
  }

  resize() {
    const c = this.canvas;
    const w = Math.max(64, Math.round(c.clientWidth * this.scale * devicePixelRatio));
    const h = Math.max(64, Math.round(c.clientHeight * this.scale * devicePixelRatio));
    if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
  }

  camera() {
    const { az, el, dist } = this.orbit;
    const cy = GEO.H * 0.48;
    const ce = Math.cos(el), se = Math.sin(el);
    const pos = [Math.sin(az) * ce * dist, cy + se * dist, Math.cos(az) * ce * dist];
    const fwd = [-pos[0], cy - pos[1], -pos[2]];
    const fl = Math.hypot(...fwd);
    for (let i = 0; i < 3; i++) fwd[i] /= fl;
    let right = [fwd[2], 0, -fwd[0]];
    const rl = Math.hypot(...right) || 1;
    for (let i = 0; i < 3; i++) right[i] /= rl;
    // up = fwd x right, not right x fwd. The other order gives a left-handed
    // basis and renders the whole lamp upside down -- which looks almost
    // plausible, because a lava lamp inverted is still a bright glow at one end
    // of a glass column.
    const up = [
      fwd[1] * right[2] - fwd[2] * right[1],
      fwd[2] * right[0] - fwd[0] * right[2],
      fwd[0] * right[1] - fwd[1] * right[0],
    ];
    return { pos, basis: [right[0], right[1], right[2], up[0], up[1], up[2], fwd[0], fwd[1], fwd[2]] };
  }

  draw(wax, lamp) {
    const gl = this.gl, v = this.volume;
    this.resize();
    const data = v.build(wax);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_3D, this.tex);
    gl.texSubImage3D(gl.TEXTURE_3D, 0, 0, 0, 0, v.nx, v.ny, v.nz,
      gl.RGBA, gl.UNSIGNED_BYTE, data);

    gl.useProgram(this.prog);
    gl.bindVertexArray(this.vao);
    gl.viewport(0, 0, this.canvas.width, this.canvas.height);
    const cam = this.camera();
    const p = this.palette;
    gl.uniform2f(this.u.uRes, this.canvas.width, this.canvas.height);
    gl.uniform3fv(this.u.uCam, cam.pos);
    gl.uniformMatrix3fv(this.u.uBasis, false, cam.basis);
    gl.uniform1f(this.u.uTan, Math.tan(0.30));
    gl.uniform3fv(this.u.uWax, p.wax);
    gl.uniform3fv(this.u.uWaxDeep, p.deep);
    gl.uniform3fv(this.u.uFluid, p.fluid);
    gl.uniform3fv(this.u.uBulb, p.bulb);
    gl.uniform1f(this.u.uSteps, this.steps);
    gl.uniform1f(this.u.uExposure, this.exposure);
    gl.uniform1f(this.u.uGlow, this.glow);
    gl.uniform1i(this.u.uView, this.view);
    gl.uniform1f(this.u.uBulbOn, this.bulbOn);
    const norm = (T) => (T - TMIN) / (TMAX - TMIN);
    gl.uniform1f(this.u.uTcross, norm(lamp ? lamp.crossover : 49));
    gl.uniform1f(this.u.uTfluidTop, norm(lamp ? lamp.column.Ttop : 40));
    gl.uniform1f(this.u.uTfluidBase, norm(lamp ? lamp.column.Tbase : 54));
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }
}
