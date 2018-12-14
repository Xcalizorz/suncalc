function getData (test_date, test_lat, test_lng) {
  test_date = new Date(test_date);
  city_information = parseFloat(test_lat) + ", " + parseFloat(test_lng) + " @ " + test_date;
  document.getElementById("city_name_header").innerHTML = city_information;

  // reformatting date for comparison API
  var day = test_date.getDate(),
      month = test_date.getMonth() + 1,
      year = test_date.getFullYear(),
      hour = test_date.getHours(),
      minutes = test_date.getMinutes();

  var new_date = [year, month, day].join('.'),
      new_time = [hour, minutes].join(':');
  
  // Comparison
  document.getElementById("getMoonPosition_link").href = "https://www.mooncalc.org/#/" + test_lat + "," + test_lng + ",16/" + new_date + "/" + new_time + "/1/0";
  document.getElementById("getPosition_link").href = "https://www.suncalc.org/#/" + test_lat + "," + test_lng + ",16/" + new_date + "/" + new_time + "/1/0";

  var moonPositionData = SunCalc.getMoonPosition(test_date, test_lat, test_lng);
  document.getElementById("getMoonPosition_Azimuth").innerHTML = moonPositionData.azimuth;
  document.getElementById("getMoonPosition_Altitude").innerHTML = moonPositionData.altitude;
  document.getElementById("getMoonPosition_Distance").innerHTML = moonPositionData.distance;
  document.getElementById("getMoonPosition_ParallacticAngle").innerHTML = moonPositionData.parallacticAngle;

  var moonIlluminationData = SunCalc.getMoonIllumination(test_date);
  document.getElementById("getMoonIllumination_fraction").innerHTML = moonIlluminationData.fraction;
  document.getElementById("getMoonIllumination_phase").innerHTML = moonIlluminationData.phase;
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
  document.getElementById("getTimes_sunrise").innerHTML = getTimesData.sunrise;
  document.getElementById("getTimes_sunriseEnd").innerHTML = getTimesData.sunriseEnd;
  document.getElementById("getTimes_sunset").innerHTML = getTimesData.sunset;
  document.getElementById("getTimes_sunsetStart").innerHTML = getTimesData.sunsetStart;

  var getMoonTimesData = SunCalc.getMoonTimes(test_date, test_lat, test_lng, false);
  document.getElementById("getMoonTimes_moonrise").innerHTML = getMoonTimesData.rise;
  document.getElementById("getMoonTimes_moonset").innerHTML = getMoonTimesData.set;

  var getPositionData = SunCalc.getPosition(test_date, test_lat, test_lng);
  document.getElementById("getPosition_azimuth").innerHTML = getPositionData.azimuth;
  document.getElementById("getPosition_altitude").innerHTML = getPositionData.altitude;

  var phaseHuntData = SunCalc.phase_hunt(date);
  document.getElementById("new_date").innerHTML = phaseHuntData.new_date;
  document.getElementById("q1_date").innerHTML = phaseHuntData.q1_date;
  document.getElementById("full_date").innerHTML = phaseHuntData.full_date;
  document.getElementById("q3_date").innerHTML = phaseHuntData.q3_date;
  document.getElementById("nextnew_date").innerHTML = phaseHuntData.nextnew_date;
  
  return false;
}