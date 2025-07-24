# Pass Prediction CLI

This simple command line tool predicts upcoming satellite passes using TLE data.
It is self contained and independent from the main `gpredict` application.

## Installation

Install dependencies using `npm install` within the `tools` directory.

```
cd tools
npm install
```

## Usage

```
node pass-predict.js --tle <file> [--min-el <degrees>] [--hours <n>] <lat> <lon> <alt>
```

- `--tle` Path to a text file containing satellite names and TLE lines in the
  standard three line format.
- `--min-el` Minimum elevation angle in degrees required for a pass to be
  reported (default `0`).
- `--hours` Number of hours to search ahead (default `24`).
- `<lat> <lon> <alt>` Observer latitude and longitude in decimal degrees and
  altitude in kilometers.

Example:

```
node pass-predict.js --tle iss.tle --min-el 30 55.0 -1.0 0.1
```
