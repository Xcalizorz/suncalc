function getData (test_date, test_lat, test_lng, test_height) {
  test_lat = parseFloat(test_lat.value);
  test_lng = parseFloat(test_lng.value);
  test_height = parseFloat(test_height.value);

  test_date = new Date(new Date(test_date.value).toUTCString());
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

  // reformatting date for comparison API
  var day = test_date.getDate(),
      month = test_date.getMonth() + 1,
      year = test_date.getFullYear(),
      hour = test_date.getHours(),
      minutes = test_date.getMinutes();

  var new_date = [year, month, day].join('.'),
      new_time = [hour, minutes].join(':');


  var sunrise = A.Coord.secondsToHMSStr(suntimes.rise);
  var sunset = A.Coord.secondsToHMSStr(suntimes.set);
  var moonrise = A.Coord.secondsToHMSStr(moontimes.rise);
  var moonset = A.Coord.secondsToHMSStr(moontimes.set)

  // recalculating with sunset time
  var sunsetArray = sunset.split(":");
  var sunsetTime = new Date(Date.UTC(year, month-1, day, sunsetArray[0], sunsetArray[1], sunsetArray[2]));

  var sunriseArray = sunrise.split(":");
  var sunriseTime = new Date(Date.UTC(year, month-1, day, sunriseArray[0], sunriseArray[1], sunriseArray[2]));
  
  var moonsetArray = moonset.split(":");
  var moonsetTime = new Date(Date.UTC(year, month-1, day, moonsetArray[0], moonsetArray[1], moonsetArray[2]));

  var moonriseArray = moonrise.split(":");
  var moonriseTime = new Date(Date.UTC(year, month-1, day, moonriseArray[0], moonriseArray[1], moonriseArray[2]));
  

  // Comparison
  document.getElementById("getMoonPosition_link").href = "https://www.mooncalc.org/#/" + test_lat + "," + test_lng + ",16/" + new_date + "/" + new_time + "/" + test_height + "/0";
  document.getElementById("getPosition_link").href = "https://www.suncalc.org/#/" + test_lat + "," + test_lng + ",16/" + new_date + "/" + new_time + "/" + test_height + "/0";

  var moonPositionData = SunCalc.getMoonPosition(test_date, test_lat, test_lng);
  document.getElementById("getMoonPosition_Azimuth").innerHTML = moontp.hz.az;
  document.getElementById("getMoonPosition_Altitude").innerHTML = moontp.hz.alt;
  document.getElementById("getMoonPosition_Distance").innerHTML = moontp.delta;
  document.getElementById("getMoonPosition_ParallacticAngle").innerHTML = moonPositionData.parallacticAngle;

  var moonIlluminationData = SunCalc.getMoonIllumination(test_date);
  document.getElementById("getMoonIllumination_fraction").innerHTML = k;
  document.getElementById("getMoonIllumination_phase").innerHTML = i;
  document.getElementById("getMoonIllumination_angle").innerHTML = moonIlluminationData.angle;

  var getTimesData = SunCalc.getTimes(test_date, test_lat, test_lng);

  document.getElementById("getTimes_solarNoon").innerHTML = getTimesData.solarNoon;
  document.getElementById("getTimes_nadir").innerHTML = getTimesData.nadir;
  document.getElementById("getTimes_dawn").innerHTML = getTimesData.dawn;
  document.getElementById("getTimes_dusk").innerHTML = getTimesData.dusk;
  document.getElementById("getTimes_goldenHour").innerHTML = getTimesData.goldenHour;
  document.getElementById("getTimes_goldenHourEnd").innerHTML = getTimesData.goldenHourEnd;
  document.getElementById("getTimes_nauticalDawn").innerHTML = getTimesData.nauticalDawn;
  document.getElementById("getTimes_nauticalDusk").innerHTML = getTimesData.nauticalDusk;
  document.getElementById("getTimes_night").innerHTML = getTimesData.night;
  document.getElementById("getTimes_nightEnd").innerHTML = getTimesData.nightEnd;
  document.getElementById("getTimes_sunrise").innerHTML = sunriseTime;
  document.getElementById("getTimes_sunriseEnd").innerHTML = getTimesData.sunriseEnd;
  document.getElementById("getTimes_sunset").innerHTML = sunsetTime;
  document.getElementById("getTimes_sunsetStart").innerHTML = getTimesData.sunsetStart;

  var getMoonTimesData = SunCalc.getMoonTimes(test_date, test_lat, test_lng, false);
  document.getElementById("getMoonTimes_moonrise").innerHTML = moonriseTime;
  document.getElementById("getMoonTimes_moonset").innerHTML = moonsetTime;

  var getPositionData = SunCalc.getPosition(test_date, test_lat, test_lng);
  document.getElementById("getPosition_azimuth").innerHTML = suntp.hz.az;
  document.getElementById("getPosition_altitude").innerHTML = suntp.hz.alt;

  var phaseHuntData = SunCalc.phase_hunt(test_date);
  document.getElementById("new_date").innerHTML = phaseHuntData.new_date;
  document.getElementById("q1_date").innerHTML = phaseHuntData.q1_date;
  document.getElementById("full_date").innerHTML = phaseHuntData.full_date;
  document.getElementById("q3_date").innerHTML = phaseHuntData.q3_date;
  document.getElementById("nextnew_date").innerHTML = phaseHuntData.nextnew_date;
  
  document.getElementById("getMoonTimes_moonset2").innerHTML =  moonsetTime;
  document.getElementById("getTimes_sunset2").innerHTML = sunsetTime;
  document.getElementById("lagTime").innerHTML = SunCalc.lagTime(moonsetTime, sunsetTime);

  var phaseData = SunCalc.phase(test_date);
  document.getElementById("p_age").innerHTML = phaseData.age;
  document.getElementById("p_angular_diameter").innerHTML = phaseData.angular_diameter;
  document.getElementById("p_moon_earth_distance").innerHTML = phaseData.distance;
  document.getElementById("p_illuminated").innerHTML = phaseData.illuminated;
  document.getElementById("p_phase").innerHTML = phaseData.phase;
  document.getElementById("p_sun_angular_diameter").innerHTML = phaseData.sun_angular_diameter;
  document.getElementById("p_sun_distance").innerHTML = phaseData.sun_distance;

  var jdo = new A.JulianDay(sunsetTime); 
  var coord = A.EclCoord.fromWgs84(test_lat, test_lng, test_height);

  var suntimes = A.Solar.times(jdo, coord);
  var moontimes = A.Moon.times(jdo, coord);

  var suntp = A.Solar.topocentricPosition(jdo, coord, true);
  var moontp = A.Moon.topocentricPosition(jdo, coord, true);
  
  var i = A.MoonIllum.phaseAngleEq2(moontp.eq, suntp.eq);
  var k = A.MoonIllum.illuminated(i);
  var j = A.MoonIllum.positionAngle(moontp.eq, suntp.eq);

  // reformatting date for comparison API
  var day = sunsetTime.getDate(),
      month = sunsetTime.getMonth() + 1,
      year = sunsetTime.getFullYear(),
      hour = sunsetTime.getHours(),
      minutes = sunsetTime.getMinutes();

  var new_date = [year, month, day].join('.'),
      new_time = [hour, minutes].join(':');

  document.getElementById("getMoonPosition_link_sunset").href = "https://www.mooncalc.org/#/" + test_lat + "," + test_lng + ",16/" + new_date + "/" + new_time + "/1/0";
  document.getElementById("getPosition_link_sunset").href = "https://www.suncalc.org/#/" + test_lat + "," + test_lng + ",16/" + new_date + "/" + new_time + "/1/0";

  document.getElementById("sunset_header").innerHTML = sunsetTime;
  
  // Moon data
  document.getElementById("getMoonPosition_azimuth_sunset").innerHTML = moontp.hz.az;
  document.getElementById("getMoonPosition_altitude_sunset").innerHTML = moontp.hz.alt;
  document.getElementById("getMoonPosition_distance_sunset").innerHTML = moontp.delta;

  // Solar data
  document.getElementById("getPosition_altitude_sunset").innerHTML = suntp.hz.alt;
  document.getElementById("getPosition_azimuth_sunset").innerHTML = suntp.hz.az;

  // Phase Calculations
  var phaseData = SunCalc.phase(sunsetTime);
  document.getElementById("p_age_sunset").innerHTML = phaseData.age;
  document.getElementById("p_angular_diameter_sunset").innerHTML = phaseData.angular_diameter;
  document.getElementById("p_moon_earth_distance_sunset").innerHTML = phaseData.distance;
  document.getElementById("p_illuminated_sunset").innerHTML = phaseData.illuminated;
  document.getElementById("p_phase_sunset").innerHTML = phaseData.phase;
  document.getElementById("p_sun_angular_diameter_sunset").innerHTML = phaseData.sun_angular_diameter;
  document.getElementById("p_sun_distance_sunset").innerHTML = phaseData.sun_distance;

  // Moon illumination
  document.getElementById("getMoonIllumination_fraction_sunset").innerHTML = k;
  document.getElementById("getMoonIllumination_phase_sunset").innerHTML = i;
  document.getElementById("getMoonIllumination_angle_sunset").innerHTML = j;

  return false;
}