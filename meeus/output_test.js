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

  return false;
}