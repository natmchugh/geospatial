--TEST--
WGS84 of Ben Nevis to OS Grid Reference sheet
--FILE--
<?php
$airyBenNevis = transform_datum(56.796556, -5.00393, GEO_WGS84, GEO_AIRY_1830);
$eastingsNorthings = coord_to_eastings_northings($airyBenNevis['lat'], $airyBenNevis['long']);
var_dump($eastingsNorthings);
var_dump(os_grid_letters($eastingsNorthings['eastings'], $eastingsNorthings['northings']));
?>
--EXPECT--
array(2) {
  ["eastings"]=>
  float(216650.01955844)
  ["northings"]=>
  float(771249.94500687)
}
string(8) "NN166712"