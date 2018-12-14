function getData (test_date, test_lat, test_lng, city_information="") {
  var getMoonPosition_link = "getMoonPosition_link",
      getPosition_link = "getPosition_link",
      getMoonPosition_Azimuth = "getMoonPosition_Azimuth",
      getMoonPosition_Altitude = "getMoonPosition_Altitude",
      getMoonPosition_Distance = "getMoonPosition_Distance",
      getMoonPosition_ParallacticAngle = "getMoonPosition_ParallacticAngle",
      getMoonIllumination_fraction = "getMoonIllumination_fraction",
      getMoonIllumination_phase = "getMoonIllumination_phase",
      getMoonIllumination_angle = "getMoonIllumination_angle",
      getTimes_solarNoon = "getTimes_solarNoon",
      getTimes_nadir = "getTimes_nadir",
      getMoonTimes_moonrise = "getMoonTimes_moonrise",
      getMoonTimes_moonset = "getMoonTimes_moonset",
      getPosition_azimuth = "getPosition_azimuth",
      getPosition_altitude = "getPosition_altitude";

  debugger;
  if (city_information==""){
    city_information = parseFloat(test_lat) + ", " + parseFloat(test_lng);
    test_date = new Date(test_date);
  }

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
  document.getElementById(getMoonPosition_link).href = "https://www.mooncalc.org/#/" + test_lat + "," + test_lng + ",16/" + new_date + "/" + new_time + "/1/0";
  document.getElementById(getPosition_link).href = "https://www.suncalc.org/#/" + test_lat + "," + test_lng + ",16/" + new_date + "/" + new_time + "/1/0";

  var moonPositionData = SunCalc.getMoonPosition(test_date, test_lat, test_lng);
  document.getElementById(getMoonPosition_Azimuth).innerHTML = moonPositionData.azimuth;
  document.getElementById(getMoonPosition_Altitude).innerHTML = moonPositionData.altitude;
  document.getElementById(getMoonPosition_Distance).innerHTML = moonPositionData.distance;
  document.getElementById(getMoonPosition_ParallacticAngle).innerHTML = moonPositionData.parallacticAngle;

  var moonIlluminationData = SunCalc.getMoonIllumination(test_date);
  document.getElementById(getMoonIllumination_fraction).innerHTML = moonIlluminationData.fraction;
  document.getElementById(getMoonIllumination_phase).innerHTML = moonIlluminationData.phase;
  document.getElementById(getMoonIllumination_angle).innerHTML = moonIlluminationData.angle;

  var getTimesData = SunCalc.getTimes(test_date, test_lat, test_lng);
  document.getElementById(getTimes_solarNoon).innerHTML = getTimesData.solarNoon;
  document.getElementById(getTimes_nadir).innerHTML = getTimesData.nadir;

  var getMoonTimesData = SunCalc.getMoonTimes(test_date, test_lat, test_lng, false);
  document.getElementById(getMoonTimes_moonrise).innerHTML = getMoonTimesData.rise;
  document.getElementById(getMoonTimes_moonset).innerHTML = getMoonTimesData.set;

  var getPositionData = SunCalc.getPosition(test_date, test_lat, test_lng);
  document.getElementById(getPosition_azimuth).innerHTML = getPositionData.azimuth;
  document.getElementById(getPosition_altitude).innerHTML = getPositionData.altitude;

  return false;
}