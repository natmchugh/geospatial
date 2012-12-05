--TEST--
WGS84 of Ben Nevis to OS Grid Reference sheet
--FILE--
<?php
// 56° 47′ 48.49″ N, 5° 0′ 21.62″ W
$lat = dms_to_decimal(56, 47, 48.49, 'N');
$long = dms_to_decimal(5, 0, 21.62, 'W');

$airyBenNevis = transform_datum($lat, $long, GEO_WGS84, GEO_AIRY_1830);
var_dump(coord_to_os_grid($airyBenNevis['lat'], $airyBenNevis['long']));
?>
--EXPECT--
NN166712