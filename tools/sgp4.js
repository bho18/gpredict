// Minimal SGP4 implementation ported from gpredict's sgp4sdp4.c
// Only supports near-earth satellites (SGP4) and necessary
// helpers for pass prediction.

const PI = Math.PI;
const TWO_PI = 2 * PI;
const DE2RA = PI / 180;
const XKE = 7.43669161e-2; // from sgp4sdp4.h
const CK2 = 5.413079e-4;
const QOMS2T = 1.880279e-9;
const S = 1.012229;
const AE = 1.0;
const XMNPDA = 1440.0;
const XKMPER = 6378.135;
const SEC_DAY = 86400.0;

function degreesToRadians(deg) { return deg * DE2RA; }
function radiansToDegrees(rad) { return rad / DE2RA; }

// Convert TLE lines to satrec structure (simplified)
function twoline2satrec(l1, l2) {
  const satrec = {};
  satrec.epochyr = parseInt(l1.substring(18, 20), 10);
  satrec.epochdays = parseFloat(l1.substring(20, 32));
  if (satrec.epochyr < 57) satrec.epochyr += 2000; else satrec.epochyr += 1900;
  const year = satrec.epochyr;
  const dayFrac = satrec.epochdays;
  const jd = jday(year, 1, 0, 0, 0, 0) + dayFrac;
  satrec.jdsatepoch = jd;
  satrec.bstar = parseFloat(l1.substring(53, 54) + '.' + l1.substring(54, 59) + 'e' + l1.substring(59, 61));
  satrec.inclo = degreesToRadians(parseFloat(l2.substring(8, 16)));
  satrec.nodeo = degreesToRadians(parseFloat(l2.substring(17, 25)));
  satrec.ecco = parseFloat('0.' + l2.substring(26, 33));
  satrec.argpo = degreesToRadians(parseFloat(l2.substring(34, 42)));
  satrec.mo = degreesToRadians(parseFloat(l2.substring(43, 51)));
  satrec.no = parseFloat(l2.substring(52, 63));
  satrec.no = satrec.no * TWO_PI / XMNPDA; // convert to rad/min
  // initialization
  sgp4init(satrec);
  return satrec;
}

// Julian day from y, m, d, h, m, s
function jday(year, month, day, hour, minute, sec) {
  return (367.0 * year - Math.floor((7 * (year + Math.floor((month + 9) / 12))) * 0.25) + Math.floor(275 * month / 9) + day + 1721013.5 + ((sec / 60 + minute) / 60 + hour) / 24);
}

function gstime(jdut1) {
  const tut1 = (jdut1 - 2451545.0) / 36525.0;
  let temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841;
  temp = (temp * Math.PI / 180) / 240.0;
  return (temp % TWO_PI);
}

function sgp4init(sat) {
  const a1 = Math.pow(XKE / sat.no, 2/3);
  const cosio = Math.cos(sat.inclo);
  const theta2 = cosio * cosio;
  const x3thm1 = 3 * theta2 - 1.0;
  const eosq = sat.ecco * sat.ecco;
  const betao2 = 1 - eosq;
  const betao = Math.sqrt(betao2);
  const del1 = 1.5 * CK2 * x3thm1 / (a1 * a1 * betao * betao2);
  const ao = a1 * (1 - del1 * (0.5 * 2/3 + del1*(1 + 134/81 * del1)));
  const delo = 1.5 * CK2 * x3thm1 / (ao*ao*betao*betao2);
  sat.no = sat.no / (1 + delo);
  sat.a = ao / (1 - delo);
  sat.perige = (sat.a*(1 - sat.ecco) - AE) * XKMPER;
  if (sat.perige < 220) sat.isSimple = true; else sat.isSimple = false;
  sat.omgdot = 0; sat.xmdot = 0; sat.xnodot = 0; // simplified
  sat.bstar = sat.bstar;
}

// Propagate using a very small subset of SGP4 (no perturbations)
function propagate(satrec, date) {
  const tsince = (date.getTime() / 60000.0) - (satrec.jdsatepoch - 2440587.5) * XMNPDA;
  const M = satrec.mo + satrec.no * tsince;
  let E = M;
  for (let i=0;i<10;i++) {
    const deltaE = (M - E + satrec.ecco*Math.sin(E)) / (1 - satrec.ecco*Math.cos(E));
    E += deltaE;
    if (Math.abs(deltaE) < 1e-8) break;
  }
  const sinE = Math.sin(E); const cosE = Math.cos(E);
  const a = satrec.a;
  const r = a*(1 - satrec.ecco*cosE);
  const v = Math.sqrt(a) * sinE / r;
  const rdot = Math.sqrt(a) * satrec.ecco * sinE / r;
  const rfdot = Math.sqrt(1 - satrec.ecco*satrec.ecco) * satrec.no / r;
  const xorb = r*(cosE - satrec.ecco);
  const yorb = r*Math.sqrt(1 - satrec.ecco*satrec.ecco)*sinE;
  const cosw = Math.cos(satrec.argpo); const sinw = Math.sin(satrec.argpo);
  const cosi = Math.cos(satrec.inclo); const sini = Math.sin(satrec.inclo);
  const cosn = Math.cos(satrec.nodeo); const sinn = Math.sin(satrec.nodeo);
  const x1 = cosw*xorb - sinw*yorb;
  const y1 = sinw*xorb + cosw*yorb;
  const x = (cosn*x1 - sinn*y1*cosi) * XKMPER;
  const y = (sinn*x1 + cosn*y1*cosi) * XKMPER;
  const z = (y1*sini) * XKMPER;
  return { position:{x,y,z} };
}

function eciToEcf(eci, gmst) {
  const x = eci.x*Math.cos(gmst) - eci.y*Math.sin(gmst);
  const y = eci.x*Math.sin(gmst) + eci.y*Math.cos(gmst);
  const z = eci.z;
  return {x,y,z};
}

function ecfToLookAngles(observer, ecf){
  const rx = ecf.x - observer.x;
  const ry = ecf.y - observer.y;
  const rz = ecf.z - observer.z;
  const topS = observer.sinlat*Math.cos(observer.lon)*rx + observer.sinlat*Math.sin(observer.lon)*ry - observer.coslat*rz;
  const topE = -Math.sin(observer.lon)*rx + Math.cos(observer.lon)*ry;
  const topZ = observer.coslat*Math.cos(observer.lon)*rx + observer.coslat*Math.sin(observer.lon)*ry + observer.sinlat*rz;
  const range = Math.sqrt(topS*topS + topE*topE + topZ*topZ);
  const elevation = Math.asin(topZ/range);
  const azimuth = Math.atan2(-topE, topS) + PI;
  return {azimuth,elevation,range};
}

module.exports = { twoline2satrec, propagate, gstime, eciToEcf, ecfToLookAngles, radiansToDegrees, degreesToRadians, jday };
