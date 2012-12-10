--TEST--
WGS84 of Sullom Voe to OS Grid Eastings and Northings
--FILE--
<?php

$airySollomVoe = transform_datum(60.459963, -1.280948, GEO_WGS84, GEO_AIRY_1830);
$eastingsNorthings = coord_to_eastings_northings($airySollomVoe['lat'], $airySollomVoe['long'], true);
var_dump(os_grid_letters($eastingsNorthings['eastings'], $eastingsNorthings['northings']));
var_dump(os_grid_numeric($eastingsNorthings['eastings'], $eastingsNorthings['northings']));
?>
--EXPECT--
string(8) "HU396753"
array(2) {
  [0]=>
  string(5) "04396"
  [1]=>
  string(5) "11750"
}