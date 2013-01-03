--TEST--
haversine() function - basic test for haversine forumla
--INI--
precision=15
--FILE--
<?php
$start = new LatLong("39.06546 N", "104.88544 W");

var_dump($start->getHaversineDistance(new LatLong("39.06546 N", "104.80 W")));
?>
--EXPECT--
float(7.38469839293155)
