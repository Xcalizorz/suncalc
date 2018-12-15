/*
 (c) 2011-2015, Vladimir Agafonkin
 SunCalc is a JavaScript library for calculating sun/moon position and light phases.
 https://github.com/mourner/suncalc
*/

(function () { 'use strict';

// shortcuts for easier to read formulas

var PI   = Math.PI,
    sin  = Math.sin,
    cos  = Math.cos,
    tan  = Math.tan,
    asin = Math.asin,
    atan = Math.atan2,
    acos = Math.acos,
    rad  = PI / 180,
    toDeg = (180 / PI);

// sun calculations are based on http://aa.quae.nl/en/reken/zonpositie.html formulas


// date/time constants and conversions

var dayMs = 1000 * 60 * 60 * 24,
    J1970 = 2440588,
    J2000 = 2451545;

function toJulian(date) { return date.valueOf() / dayMs - 0.5 + J1970; }
function fromJulian(j)  { return new Date((j + 0.5 - J1970) * dayMs); }
function toDays(date)   { return toJulian(date) - J2000; }


// general calculations for position
function toDegree(radians) {
  return radians * 180 / Math.PI;
}

var e = rad * 23.4397; // obliquity of the Earth

function rightAscension(l, b) { return atan(sin(l) * cos(e) - tan(b) * sin(e), cos(l)); }
function declination(l, b)    { return asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l)); }

function azimuth(H, phi, dec)  { return 180 + toDegree(atan(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi))); }
function altitude(H, phi, dec) { return toDegree(asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H))); }

function siderealTime(d, lw) { return rad * (280.16 + 360.9856235 * d) - lw; }

function astroRefraction(h) {
    if (h < 0) // the following formula works for positive altitudes only.
        h = 0; // if h = -0.08901179 a div/0 would occur.

    // formula 16.4 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
    // 1.02 / tan(h + 10.26 / (h + 5.10)) h in degrees, result in arc minutes -> converted to rad:
    return 0.0002967 / Math.tan(h + 0.00312536 / (h + 0.08901179));
}

// general sun calculations

function solarMeanAnomaly(d) { return rad * (357.5291 + 0.98560028 * d); }

function eclipticLongitude(M) {

    var C = rad * (1.9148 * sin(M) + 0.02 * sin(2 * M) + 0.0003 * sin(3 * M)), // equation of center
        P = rad * 102.9372; // perihelion of the Earth

    return M + C + P + PI;
}

function sunCoords(d) {

    var M = solarMeanAnomaly(d),
        L = eclipticLongitude(M);

    return {
        dec: declination(L, 0),
        ra: rightAscension(L, 0)
    };
}


var SunCalc = {};


// calculates sun position for a given date and latitude/longitude

SunCalc.getPosition = function (date, lat, lng) {

    var lw  = rad * -lng,
        phi = rad * lat,
        d   = toDays(date),

        c  = sunCoords(d),
        H  = siderealTime(d, lw) - c.ra;

    return {
        azimuth: azimuth(H, phi, c.dec),
        altitude: altitude(H, phi, c.dec)
    };
};


// sun times configuration (angle, morning name, evening name)

var times = SunCalc.times = [
    [-0.833, 'sunrise',       'sunset'      ],
    [  -0.3, 'sunriseEnd',    'sunsetStart' ],
    [    -6, 'dawn',          'dusk'        ],
    [   -12, 'nauticalDawn',  'nauticalDusk'],
    [   -18, 'nightEnd',      'night'       ],
    [     6, 'goldenHourEnd', 'goldenHour'  ]
];

// adds a custom time to the times config

SunCalc.addTime = function (angle, riseName, setName) {
    times.push([angle, riseName, setName]);
};


// calculations for sun times

var J0 = 0.0009;

function julianCycle(d, lw) { return Math.round(d - J0 - lw / (2 * PI)); }

function approxTransit(Ht, lw, n) { return J0 + (Ht + lw) / (2 * PI) + n; }
function solarTransitJ(ds, M, L)  { return J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L); }

function hourAngle(h, phi, d) { return acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d))); }

// returns set time for the given sun altitude
function getSetJ(h, lw, phi, dec, n, M, L) {

    var w = hourAngle(h, phi, dec),
        a = approxTransit(w, lw, n);
    return solarTransitJ(a, M, L);
}


// calculates sun times for a given date and latitude/longitude

SunCalc.getTimes = function (date, lat, lng) {

    var lw = rad * -lng,
        phi = rad * lat,

        d = toDays(date),
        n = julianCycle(d, lw),
        ds = approxTransit(0, lw, n),

        M = solarMeanAnomaly(ds),
        L = eclipticLongitude(M),
        dec = declination(L, 0),

        Jnoon = solarTransitJ(ds, M, L),

        i, len, time, Jset, Jrise;


    var result = {
        solarNoon: fromJulian(Jnoon),
        nadir: fromJulian(Jnoon - 0.5)
    };

    for (i = 0, len = times.length; i < len; i += 1) {
        time = times[i];

        Jset = getSetJ(time[0] * rad, lw, phi, dec, n, M, L);
        Jrise = Jnoon - (Jset - Jnoon);

        result[time[1]] = fromJulian(Jrise);
        result[time[2]] = fromJulian(Jset);
    }

    return result;
};


// moon calculations, based on http://aa.quae.nl/en/reken/hemelpositie.html formulas

function moonCoords(d) { // geocentric ecliptic coordinates of the moon

    var L = rad * (218.316 + 13.176396 * d), // ecliptic longitude
        M = rad * (134.963 + 13.064993 * d), // mean anomaly
        F = rad * (93.272 + 13.229350 * d),  // mean distance

        l  = L + rad * 6.289 * sin(M), // longitude
        b  = rad * 5.128 * sin(F),     // latitude
        dt = 385001 - 20905 * cos(M);  // distance to the moon in km

    return {
        ra: rightAscension(l, b),
        dec: declination(l, b),
        dist: dt
    };
}

SunCalc.getMoonPosition = function (date, lat, lng) {

    var lw  = rad * -lng,
        phi = rad * lat,
        d   = toDays(date),

        c = moonCoords(d),
        H = siderealTime(d, lw) - c.ra,
        h = altitude(H, phi, c.dec),
        // formula 14.1 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
        pa = atan(sin(H), tan(phi) * cos(c.dec) - sin(c.dec) * cos(H));

    h = h + astroRefraction(h); // altitude correction for refraction

    return {
        azimuth: azimuth(H, phi, c.dec),
        altitude: h,
        distance: c.dist,
        parallacticAngle: pa
    };
};


// calculations for illumination parameters of the moon,
// based on http://idlastro.gsfc.nasa.gov/ftp/pro/astro/mphase.pro formulas and
// Chapter 48 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.

SunCalc.getMoonIllumination = function (date) {

    var d = toDays(date || new Date()),
        s = sunCoords(d),
        m = moonCoords(d),

        // according to https://www.space.com/17081-how-far-is-earth-from-the-sun.html
        sdist = 149597870, // distance from Earth to Sun in km

        phi = acos(sin(s.dec) * sin(m.dec) + cos(s.dec) * cos(m.dec) * cos(s.ra - m.ra)),
        inc = atan(sdist * sin(phi), m.dist - sdist * cos(phi)),
        angle = atan(cos(s.dec) * sin(s.ra - m.ra), sin(s.dec) * cos(m.dec) -
                cos(s.dec) * sin(m.dec) * cos(s.ra - m.ra)),
        phase_original = 0.5 + 0.5 * inc * (angle < 0 ? -1 : 1) / Math.PI,
        phase = parseFloat(phase_original).toFixed(8) + ", ";

        switch(true) {
          case phase_original == 0:
            phase += "New Moon";
            break;
          case phase_original < 0.25:
            phase += "Waxing Gibbous";
            break;
          case phase_original == 0.25:
            phase += "First Quarter";
            break;
          case phase_original < 0.50:
            phase += "Waxing Gibbous";
            break;
          case phase_original == 0.50:
            phase += "Full Moon";
            break;
          case phase_original < 0.75:
            phase += "Waning Gibbous";
            break;
          case phase_original == 0.75:
            phase += "Last Quarter";
            break;
          case phase_original <= 1.00:
            phase += "Waning Crescent";
            break;
          default:
            phase = "Unknown";
        }

    return {
        fraction: parseFloat((((1 + cos(inc)) / 2) * 100).toFixed(4)),
        phase: phase,
        angle: toDegree(angle)
    };
};


function hoursLater(date, h) {
    return new Date(date.valueOf() + h * dayMs / 24);
}

// calculations for moon rise/set times are based on http://www.stargazing.net/kepler/moonrise.html article

SunCalc.getMoonTimes = function (date, lat, lng, inUTC) {
    var t = new Date(date);
    if (inUTC) t.setUTCHours(0, 0, 0, 0);
    else t.setHours(0, 0, 0, 0);

    var hc = 0.133 * rad,
        h0 = SunCalc.getMoonPosition(t, lat, lng).altitude - hc,
        h1, h2, rise, set, a, b, xe, ye, d, roots, x1, x2, dx;

    // go in 2-hour chunks, each time seeing if a 3-point quadratic curve crosses zero (which means rise or set)
    for (var i = 1; i <= 24; i += 2) {
        h1 = SunCalc.getMoonPosition(hoursLater(t, i), lat, lng).altitude - hc;
        h2 = SunCalc.getMoonPosition(hoursLater(t, i + 1), lat, lng).altitude - hc;

        a = (h0 + h2) / 2 - h1;
        b = (h2 - h0) / 2;
        xe = -b / (2 * a);
        ye = (a * xe + b) * xe + h1;
        d = b * b - 4 * a * h1;
        roots = 0;

        if (d >= 0) {
            dx = Math.sqrt(d) / (Math.abs(a) * 2);
            x1 = xe - dx;
            x2 = xe + dx;
            if (Math.abs(x1) <= 1) roots++;
            if (Math.abs(x2) <= 1) roots++;
            if (x1 < -1) x1 = x2;
        }

        if (roots === 1) {
            if (h0 < 0) rise = i + x1;
            else set = i + x1;

        } else if (roots === 2) {
            rise = i + (ye < 0 ? x2 : x1);
            set = i + (ye < 0 ? x1 : x2);
        }

        if (rise && set) break;

        h0 = h2;
    }

    var result = {};

    if (rise) result.rise = hoursLater(t, rise);
    if (set) result.set = hoursLater(t, set);

    if (!rise && !set) result[ye > 0 ? 'alwaysUp' : 'alwaysDown'] = true;

    return result;
};

/** Julian Calculations */
function fromDate (date) {
  return date.getTime() / 86400000 + 2440587.5
}

function toDate (julian) {
  return new Date((julian - 2440587.5) * 86400000)
}


// Phases of the moon & precision
var NEW = 0
var FIRST = 1
var FULL = 2
var LAST = 3
var PHASE_MASK = 3

// Astronomical Constants
// JDN stands for Julian Day Number
// Angles here are in degrees
// 1980 January 0.0 in JDN
// XXX: DateTime(1980).jdn yields 2444239.5 -- which one is right?
// XXX: even though 2444239.5 is correct for the 1 Jan 1980, 2444238.5 gives
// better accuracy results... possibly somebody chose all of the below
// constants based on the wrong epoch?
const EPOCH = 2444238.5

// Ecliptic longitude of the Sun at epoch 1980.0
const ECLIPTIC_LONGITUDE_EPOCH = 278.833540

// Ecliptic longitude of the Sun at perigee
const ECLIPTIC_LONGITUDE_PERIGEE = 282.596403

// Eccentricity of Earth's orbit
const ECCENTRICITY = 0.016718

// Semi-major axis of Earth's orbit, in kilometers
const SUN_SMAXIS = 1.49585e8

// Sun's angular size, in degrees, at semi-major axis distance
const SUN_ANGULAR_SIZE_SMAXIS = 0.533128

// Elements of the Moon's orbit, epoch 1980.0
// Moon's mean longitude at the epoch
const MOON_MEAN_LONGITUDE_EPOCH = 64.975464

// Mean longitude of the perigee at the epoch
const MOON_MEAN_PERIGEE_EPOCH = 349.383063

// Eccentricity of the Moon's orbit
const MOON_ECCENTRICITY = 0.054900

// Semi-major axis of the Moon's orbit, in kilometers
const MOON_SMAXIS = 384401.0

// MOON_SMAXIS premultiplied by the angular size of the Moon from the Earth
const MOON_ANGULAR_SIZE_SMAXIS = MOON_SMAXIS * 0.5181

// Synodic month (new Moon to new Moon), in days
const SYNODIC_MONTH = 29.53058868


// sin cos functions
function dsin (d) {
  return Math.sin(torad(d))
}

function dcos (d) {
  return Math.cos(torad(d))
}

/**
 * Calculates time of the mean new Moon for a given base date.
 * This argument K to this function is the precomputed synodic month
 * index, given by:
 *   K = (year - 1900) * 12.3685
 * where year is expressed as a year and fractional year.
 * @param  {Date} sdate   Start date
 * @param  {[type]} k     [description]
 * @return {[type]}       [description]
 */
function meanphase (sdate, k) {
  // Time in Julian centuries from 1900 January 12 noon UTC
  var delta_t = (sdate - -2208945600000.0) / 86400000.0;
  var t = delta_t / 36525;
  return 2415020.75933 +
    SYNODIC_MONTH * k +
    (0.0001178 - 0.000000155 * t) * t * t +
    0.00033 * dsin(166.56 + (132.87 - 0.009173 * t) * t);
}

/**
 * Given a K value used to determine the mean phase of the new moon, and a
 * phase selector (0, 1, 2, 3), obtain the true, corrected phase time.
 * @param  {[type]} k      [description]
 * @param  {[type]} tphase [description]
 * @return {[type]}        [description]
 */
function truephase (k, tphase) {
  // restrict tphase to (0, 1, 2, 3)
  tphase = tphase & PHASE_MASK;

  // add phase to new moon time
  k = k + 0.25 * tphase;

  // Time in Julian centuries from 1900 January 0.5
  var t = (1.0 / 1236.85) * k;

  // Mean time of phase
  var pt = 2415020.75933 +
    SYNODIC_MONTH * k +
    (0.0001178 - 0.000000155 * t) * t * t +
    0.00033 * dsin(166.56 + (132.87 - 0.009173 * t) * t);

  // Sun's mean anomaly
  var m = 359.2242 + 29.10535608 * k - (0.0000333 - 0.00000347 * t) * t * t;

  // Moon's mean anomaly
  var mprime = 306.0253 + 385.81691806 * k + (0.0107306 + 0.00001236 * t) * t * t;

  // Moon's argument of latitude
  var f = 21.2964 + 390.67050646 * k - (0.0016528 - 0.00000239 * t) * t * t;

  // use different correction equations depending on the phase being sought
  switch (tphase) {
    // new and full moon use one correction
    case NEW:
    case FULL:
      pt += (0.1734 - 0.000393 * t) * dsin(m) +
        0.0021 * dsin(2 * m) -
        0.4068 * dsin(mprime) +
        0.0161 * dsin(2 * mprime) -
        0.0004 * dsin(3 * mprime) +
        0.0104 * dsin(2 * f) -
        0.0051 * dsin(m + mprime) -
        0.0074 * dsin(m - mprime) +
        0.0004 * dsin(2 * f + m) -
        0.0004 * dsin(2 * f - m) -
        0.0006 * dsin(2 * f + mprime) +
        0.0010 * dsin(2 * f - mprime) +
        0.0005 * dsin(m + 2 * mprime);
      break;

    // first and last quarter moon use a different correction
    case FIRST:
    case LAST:
      pt += (0.1721 - 0.0004 * t) * dsin(m) +
        0.0021 * dsin(2 * m) -
        0.6280 * dsin(mprime) +
        0.0089 * dsin(2 * mprime) -
        0.0004 * dsin(3 * mprime) +
        0.0079 * dsin(2 * f) -
        0.0119 * dsin(m + mprime) -
        0.0047 * dsin(m - mprime) +
        0.0003 * dsin(2 * f + m) -
        0.0004 * dsin(2 * f - m) -
        0.0006 * dsin(2 * f + mprime) +
        0.0021 * dsin(2 * f - mprime) +
        0.0003 * dsin(m + 2 * mprime) +
        0.0004 * dsin(m - 2 * mprime) -
        0.0003 * dsin(2 * m + mprime);

      // the sign of the last term depends on whether we're looking for a first
      // or last quarter moon!
      var sign = (tphase < FULL) ? +1 : -1;
      pt += sign * (0.0028 - 0.0004 * dcos(m) + 0.0003 * dcos(mprime));

      break;
  }

  return toDate(pt);
}

/**
 * Find time of phases of the moon which surround the current date.
 * Five phases are found, starting and ending with the new moons
 * which bound the current lunation.
 * @param  {Date} sdate Date to start hunting from (defaults to current date)
 * @return {Object}     Object containing recent past and future phases
 */
SunCalc.phase_hunt = function (sdate) {
  if (!sdate) {
    sdate = new Date();
  }

  var adate = new Date(sdate.getTime() - (45 * 86400000)); // 45 days prior
  var k1 = Math.floor(12.3685 * (adate.getFullYear() + (1.0 / 12.0) * adate.getMonth() - 1900));
  var nt1 = meanphase(adate.getTime(), k1);

  sdate = fromDate(sdate);
  adate = nt1 + SYNODIC_MONTH;
  var k2 = k1 + 1;
  var nt2 = meanphase(adate, k2);
  while (nt1 > sdate || sdate >= nt2) {
    adate += SYNODIC_MONTH;
    k1++;
    k2++;
    nt1 = nt2;
    nt2 = meanphase(adate, k2);
  }

  return {
    new_date: truephase(k1, NEW),
    q1_date: truephase(k1, FIRST),
    full_date: truephase(k1, FULL),
    q3_date: truephase(k1, LAST),
    nextnew_date: truephase(k2, NEW)
  }
}


/**
 * Finds the phase information for specific date.
 * @param  {Date} phase_date Date to get phase information of.
 * @return {Object}          Phase data
 */
SunCalc.phase = function (phase_date) {
  if (!phase_date) {
    phase_date = new Date()
  }
  phase_date = fromDate(phase_date)

  const day = phase_date - EPOCH

  // calculate sun position
  const sun_mean_anomaly =
    (360.0 / 365.2422) * day +
    (ECLIPTIC_LONGITUDE_EPOCH - ECLIPTIC_LONGITUDE_PERIGEE)
  const sun_true_anomaly =
    2 * toDegree(Math.atan(
      Math.sqrt((1.0 + ECCENTRICITY) / (1.0 - ECCENTRICITY)) *
      Math.tan(0.5 * kepler(sun_mean_anomaly, ECCENTRICITY))
    ))
  const sun_ecliptic_longitude =
    ECLIPTIC_LONGITUDE_PERIGEE + sun_true_anomaly
  const sun_orbital_distance_factor =
    (1 + ECCENTRICITY * dcos(sun_true_anomaly)) /
    (1 - ECCENTRICITY * ECCENTRICITY)

  // calculate moon position
  const moon_mean_longitude =
    MOON_MEAN_LONGITUDE_EPOCH + 13.1763966 * day
  const moon_mean_anomaly =
    moon_mean_longitude - 0.1114041 * day - MOON_MEAN_PERIGEE_EPOCH
  const moon_evection =
    1.2739 * dsin(
      2 * (moon_mean_longitude - sun_ecliptic_longitude) - moon_mean_anomaly
    )
  const moon_annual_equation =
    0.1858 * dsin(sun_mean_anomaly)
  // XXX: what is the proper name for this value?
  const moon_mp =
    moon_mean_anomaly +
    moon_evection -
    moon_annual_equation -
    0.37 * dsin(sun_mean_anomaly)
  const moon_equation_center_correction =
    6.2886 * dsin(moon_mp)
  const moon_corrected_longitude =
    moon_mean_longitude +
    moon_evection +
    moon_equation_center_correction -
    moon_annual_equation +
    0.214 * dsin(2.0 * moon_mp)
  const moon_age =
    fixangle(
      moon_corrected_longitude -
      sun_ecliptic_longitude +
      0.6583 * dsin(
        2 * (moon_corrected_longitude - sun_ecliptic_longitude)
      )
    )
  const moon_distance =
    (MOON_SMAXIS * (1.0 - MOON_ECCENTRICITY * MOON_ECCENTRICITY)) /
    (1.0 + MOON_ECCENTRICITY * dcos(moon_mp + moon_equation_center_correction))

  return {
    phase: (1.0 / 360.0) * moon_age,
    illuminated: 0.5 * (1.0 - dcos(moon_age)),
    age: (SYNODIC_MONTH / 360.0) * moon_age,
    distance: moon_distance,
    angular_diameter: MOON_ANGULAR_SIZE_SMAXIS / moon_distance,
    sun_distance: SUN_SMAXIS / sun_orbital_distance_factor,
    sun_angular_diameter: SUN_ANGULAR_SIZE_SMAXIS * sun_orbital_distance_factor
  }
}

/**
 * Returns a HH:MM:SS formatted string from 
 * a given milliseconds input
 * 
 * @param {Int} msec_num Milliseconds 
 */
function toHHMMSS (msec_num) {
  var sec_num = msec_num / 1000;

  var days = Math.floor(sec_num / (3600*24))
  var hours   = Math.floor((sec_num - (days * 3600 * 24)) / 3600);
  var minutes = Math.floor((sec_num - (days * 3600 * 24) - (hours * 3600)) / 60);
  var seconds = sec_num - (days * 3600 * 24) - (hours * 3600) - (minutes * 60);
  
  if (days < 0) {
    var daysHourDifference = days*24 + hours;
    days -= Math.floor(daysHourDifference / 24);
    hours = Math.floor((daysHourDifference - (days * 24)));
  }

  seconds = seconds.toFixed(2);

  if (0 <= hours && hours < 10) {hours   = "0"+hours;}
  if (minutes < 10) {minutes = "0"+minutes;}
  if (seconds < 10) {seconds = "0"+seconds;}
  return days+' days, '+hours+':'+minutes+':'+seconds;
}

/**
 * The time between moonset and sunset
 * @param moonset Time of moonset 
 *    Calculate via Sunrise.getMoonTimes(..).moonset
 * @param sunset Time of sunset
 *    Calculate via Sunrise.getTimes(..).sunset
 */
SunCalc.lagTime = function (moonset, sunset) {
  return toHHMMSS(moonset-sunset);
}

// export as Node module / AMD module / browser variable
if (typeof exports === 'object' && typeof module !== 'undefined') module.exports = SunCalc;
else if (typeof define === 'function' && define.amd) define(SunCalc);
else window.SunCalc = SunCalc;

}());
