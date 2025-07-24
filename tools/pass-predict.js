#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const sgp4 = require('./sgp4');
const yargs = require('yargs');
const { hideBin } = require('yargs/helpers');

const argv = yargs(hideBin(process.argv))
  .option('tle', {
    alias: 't',
    describe: 'Path to TLE file containing pairs of lines per satellite',
    demandOption: true,
    type: 'string'
  })
  .option('min-el', {
    alias: 'm',
    describe: 'Minimum elevation in degrees for passes',
    type: 'number',
    default: 0
  })
  .option('hours', {
    alias: 'h',
    describe: 'Number of hours to search ahead',
    type: 'number',
    default: 24
  })
  .demandCommand(3, 'Please provide latitude, longitude, and altitude')
  .usage('Usage: $0 [options] <lat> <lon> <alt>')
  .help()
  .argv;

const lat = parseFloat(argv._[0]);
const lon = parseFloat(argv._[1]);
const alt = parseFloat(argv._[2]);

if (isNaN(lat) || isNaN(lon) || isNaN(alt)) {
  console.error('Invalid latitude, longitude, or altitude');
  process.exit(1);
}

function readTLEs(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').split(/\r?\n/).filter(l => l.trim());
  const tles = [];
  for (let i = 0; i < lines.length; i += 3) {
    const name = lines[i].trim();
    const l1 = lines[i + 1];
    const l2 = lines[i + 2];
    if (l1 && l2) {
      tles.push({ name, tle1: l1.trim(), tle2: l2.trim() });
    }
  }
  return tles;
}

function degrees(radians) {
  return radians * 180 / Math.PI;
}

function formatDate(date) {
  return date.toISOString().replace(/T/, ' ').replace(/\..+/, '');
}

function predictPasses(tle, obs, hours, minEl) {
  const satrec = sgp4.twoline2satrec(tle.tle1, tle.tle2);
  const start = new Date();
  const end = new Date(start.getTime() + hours * 3600 * 1000);
  let time = new Date(start);
  const stepMs = 10000; // 10 seconds
  let above = false;
  let aos = null;
  let aosAz = 0;
  let maxEl = -90;
  let maxElAz = 0;
  const passes = [];

  while (time <= end) {
    const eci = sgp4.propagate(satrec, time);
    if (!eci.position) {
      time = new Date(time.getTime() + stepMs);
      continue;
    }
    const gmst = sgp4.gstime(sgp4.jday(time.getUTCFullYear(), time.getUTCMonth()+1, time.getUTCDate(), time.getUTCHours(), time.getUTCMinutes(), time.getUTCSeconds()));
    const ecf = sgp4.eciToEcf(eci.position, gmst);
    const look = sgp4.ecfToLookAngles(obs, ecf);
    const elev = degrees(look.elevation);
    const az = (degrees(look.azimuth) + 360) % 360;

    if (elev > maxEl) {
      maxEl = elev;
      maxElAz = az;
    }

    if (!above && elev >= 0) {
      // AOS
      above = true;
      aos = new Date(time);
      aosAz = az;
      maxEl = elev;
      maxElAz = az;
    } else if (above && elev < 0) {
      // LOS
      const los = new Date(time);
      if (maxEl >= minEl) {
        passes.push({ aos, aosAz, los, losAz: az, maxEl, maxElAz });
      }
      above = false;
      maxEl = -90;
    }

    time = new Date(time.getTime() + stepMs);
  }

  return passes;
}

const observer = {
  longitude: sgp4.degreesToRadians(lon),
  latitude: sgp4.degreesToRadians(lat),
  height: alt
};

const tles = readTLEs(path.resolve(argv.tle));
if (tles.length === 0) {
  console.error('No TLEs found');
  process.exit(1);
}

for (const tle of tles) {
  const passes = predictPasses(tle, observer, argv.hours, argv['min-el']);
  if (passes.length === 0) {
    console.log(`No passes found for ${tle.name}`);
    continue;
  }
  console.log(`\nPasses for ${tle.name}`);
  console.log('AOS (UTC)            Az   LOS (UTC)            Az   MaxEl');
  for (const p of passes) {
    console.log(`${formatDate(p.aos)} ${p.aosAz.toFixed(1).padStart(5)}  ${formatDate(p.los)} ${p.losAz.toFixed(1).padStart(5)}  ${p.maxEl.toFixed(1).padStart(6)}`);
  }
}
