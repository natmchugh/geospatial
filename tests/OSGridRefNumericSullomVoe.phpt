--TEST--
WGS84 of Sullom Voe to OS Grid Eastings and Northings
--FILE--
<?php
// 56° 47′ 48.49″ N, 5° 0′ 21.62″ W
$lat = dms_to_decimal(56, 47, 48.49, 'N');
$long = dms_to_decimal(5, 0, 21.62, 'W');
$airySollomVoe = transform_datum($lat, $long, GEO_WGS84, GEO_AIRY_1830);
var_dump($airySollomVoe);
var_dump(coord_to_os_grid($airySollomVoe['lat'], $airySollomVoe['long'], true));
?>
--EXPECT--
array(2) {
  ["eastings"]=>
  int(04396)
  ["northings"]=>
  int(11753)
}