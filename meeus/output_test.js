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
var EPOCH = 2444238.5

// Ecliptic longitude of the Sun at epoch 1980.0
var ECLIPTIC_LONGITUDE_EPOCH = 278.833540

// Ecliptic longitude of the Sun at perigee
var ECLIPTIC_LONGITUDE_PERIGEE = 282.596403

// Eccentricity of Earth's orbit
var ECCENTRICITY = 0.016718

// Semi-major axis of Earth's orbit, in kilometers
var SUN_SMAXIS = 1.49585e8

// Sun's angular size, in degrees, at semi-major axis distance
var SUN_ANGULAR_SIZE_SMAXIS = 0.533128

// Elements of the Moon's orbit, epoch 1980.0
// Moon's mean longitude at the epoch
var MOON_MEAN_LONGITUDE_EPOCH = 64.975464

// Mean longitude of the perigee at the epoch
var MOON_MEAN_PERIGEE_EPOCH = 349.383063

// Eccentricity of the Moon's orbit
var MOON_ECCENTRICITY = 0.054900

// Semi-major axis of the Moon's orbit, in kilometers
var MOON_SMAXIS = 384401.0

// MOON_SMAXIS premultiplied by the angular size of the Moon from the Earth
var MOON_ANGULAR_SIZE_SMAXIS = MOON_SMAXIS * 0.5181

// Synodic month (new Moon to new Moon), in days
var SYNODIC_MONTH = 29.53058868


// sin cos functions
function dsin (d) {
  return Math.sin(torad(d))
}

function dcos (d) {
  return Math.cos(torad(d))
}

/**
 * Solve the equation of Kepler.
 */
function kepler (m, ecc) {
  var epsilon = 1e-6;

  m = torad(m);
  var e = m;
  while (1) {
    var delta = e - ecc * Math.sin(e) - m;
    e -= delta / (1.0 - ecc * Math.cos(e));

    if (Math.abs(delta) <= epsilon) {
      break;
    }
  }

  return e;
}

function fixangle (a) {
  return a - 360.0 * Math.floor(a / 360.0);
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
phase_hunt = function (sdate) {
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
phase = function (phase_date) {
  if (!phase_date) {
    phase_date = new Date();
  }
  phase_date = fromDate(phase_date);

  var day = phase_date - EPOCH;

  // calculate sun position
  var sun_mean_anomaly =
    (360.0 / 365.2422) * day +
    (ECLIPTIC_LONGITUDE_EPOCH - ECLIPTIC_LONGITUDE_PERIGEE);
  var sun_true_anomaly =
    2 * todeg(
      Math.atan(Math.sqrt((1.0 + ECCENTRICITY) / (1.0 - ECCENTRICITY)) *
      Math.tan(0.5 * kepler(sun_mean_anomaly, ECCENTRICITY)))
    );
  var sun_ecliptic_longitude =
    ECLIPTIC_LONGITUDE_PERIGEE + sun_true_anomaly;
  var sun_orbital_distance_factor =
    (1 + ECCENTRICITY * dcos(sun_true_anomaly)) /
    (1 - ECCENTRICITY * ECCENTRICITY);

  // calculate moon position
  var moon_mean_longitude =
    MOON_MEAN_LONGITUDE_EPOCH + 13.1763966 * day;
  var moon_mean_anomaly =
    moon_mean_longitude - 0.1114041 * day - MOON_MEAN_PERIGEE_EPOCH;
  var moon_evection =
    1.2739 * dsin(
      2 * (moon_mean_longitude - sun_ecliptic_longitude) - moon_mean_anomaly
    );
  var moon_annual_equation =
    0.1858 * dsin(sun_mean_anomaly);
  // XXX: what is the proper name for this value?
  var moon_mp =
    moon_mean_anomaly +
    moon_evection -
    moon_annual_equation -
    0.37 * dsin(sun_mean_anomaly);
  var moon_equation_center_correction =
    6.2886 * dsin(moon_mp);
  var moon_corrected_longitude =
    moon_mean_longitude +
    moon_evection +
    moon_equation_center_correction -
    moon_annual_equation +
    0.214 * dsin(2.0 * moon_mp);
  var moon_age =
    fixangle(
      moon_corrected_longitude -
      sun_ecliptic_longitude +
      0.6583 * dsin(
        2 * (moon_corrected_longitude - sun_ecliptic_longitude)
      )
    );
  var moon_distance =
    (MOON_SMAXIS * (1.0 - MOON_ECCENTRICITY * MOON_ECCENTRICITY)) /
    (1.0 + MOON_ECCENTRICITY * dcos(moon_mp + moon_equation_center_correction));

  return {
    phase: (1.0 / 360.0) * moon_age,
    illuminated: 0.5 * (1.0 - dcos(moon_age)),
    age: (SYNODIC_MONTH / 360.0) * moon_age,
    distance: moon_distance,
    angular_diameter: MOON_ANGULAR_SIZE_SMAXIS / moon_distance,
    sun_distance: SUN_SMAXIS / sun_orbital_distance_factor,
    sun_angular_diameter: SUN_ANGULAR_SIZE_SMAXIS * sun_orbital_distance_factor
  };
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
lagTime = function (moonset, sunset) {
  return toHHMMSS(moonset-sunset);
}




function getData (test_date, test_lat, test_lng, test_height) {
  test_lat = parseFloat(test_lat.value);
  test_lng = parseFloat(test_lng.value);
  test_height = parseFloat(test_height.value);

  test_date = new Date(test_date.value);
  city_information = parseFloat(test_lat) + ", " + parseFloat(test_lng) + " @ " + test_date;
  document.getElementById("city_name_header").innerHTML = city_information;

  var jdo = new A.JulianDay(test_date); 
  var coord = A.EclCoord.fromWgs84(test_lat, test_lng, test_height);

  var suntimes = A.Solar.times(jdo, coord);
  var moontimes = A.Moon.times(jdo, coord);

  var suntp = A.Solar.topocentricPosition(jdo, coord, true);
  var moontp = A.Moon.topocentricPosition(jdo, coord, true);
  
  var i = A.MoonIllum.phaseAngleEq2(moontp.eq, suntp.eq);
  var k = A.MoonIllum.illuminated(i);
  
  $('#suntimes').text("rise: " + A.Coord.secondsToHMSStr(suntimes.rise) + 
            ", transit: " + A.Coord.secondsToHMSStr(suntimes.transit) + 
            ", set: " +  A.Coord.secondsToHMSStr(suntimes.set));
            
  $('#moontimes').text("rise: " + A.Coord.secondsToHMSStr(moontimes.rise) + 
            ", transit: " + A.Coord.secondsToHMSStr(moontimes.transit) + 
            ", set: " +  A.Coord.secondsToHMSStr(moontimes.set));

  $('#date').text(test_date.toString() + ", jd: " + jdo.jd);
  $('#sunpos').text(suntp.hz.toString());
  $('#moonpos').text(moontp.hz.toString() + ", dist: " + moontp.delta);		
  
  $('#moonillum').text("phase: " + i + ", illuminated: " + k);

  var phaseHuntData = phase_hunt(test_date);
  $('#moondata').text("Most Recent New Moon : " + phaseHuntData.new_date +
                      "First Quarter: " + phaseHuntData.q1_date +
                      "Full Moon: " + phaseHuntData.full_date +
                      "3rd Quarter: " + phaseHuntData.q3_date +
                      "Next New Moon: " + phaseHuntData.nextnew_date
  );

  return false;
}