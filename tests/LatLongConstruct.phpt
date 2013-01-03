--TEST--
Test LatLong Constructor
--FILE--
<?php

$latLong = new LatLong("52.5 N", "1.2 W");
var_dump($latLong);

$latLong2 = new LatLong("52.5", "-1.2");
var_dump($latLong2);

$latLong3 = new LatLong("-52.5 S", "1.2 W");
var_dump($latLong3);

$latLong4 = new LatLong('-52° 30.0" S', "1.2 W");
var_dump($latLong4);
--EXPECT--
object(LatLong)#1 (3) {
  ["latitude"]=>
  float(52.5)
  ["longitude"]=>
  float(-1.2)
  ["reference_ellipsoid"]=>
  int(1)
}
object(LatLong)#2 (3) {
  ["latitude"]=>
  float(52.5)
  ["longitude"]=>
  float(-1.2)
  ["reference_ellipsoid"]=>
  int(1)
}
object(LatLong)#3 (3) {
  ["latitude"]=>
  float(52.5)
  ["longitude"]=>
  float(-1.2)
  ["reference_ellipsoid"]=>
  int(1)
}
object(LatLong)#4 (3) {
  ["latitude"]=>
  float(-52)
  ["longitude"]=>
  float(-1.2)
  ["reference_ellipsoid"]=>
  int(1)
}