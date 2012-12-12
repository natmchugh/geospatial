--TEST--
Great Yarmouth Example From OS docs from Eastings and Northings  to lat long
--FILE--
<?php

$latLong = eastings_northings_to_coords(651409.903, 313177.270, true);
var_dump($latLong);
var_dump(decimal_to_dms($latLong['latitude']), decimal_to_dms($latLong['longitude']));
?>
--EXPECT--
array(14) {
  ["latitude"]=>
  float(52.657570301519)
  ["longitude"]=>
  float(1.7179215809789)
  ["phiDerivative"]=>
  float(0.92006620953394)
  ["M"]=>
  float(413177.26996734)
  ["nu"]=>
  float(6388523.3414779)
  ["rho"]=>
  float(6372819.3093685)
  ["eta2"]=>
  float(0.0024642205195358)
  ["VII"]=>
  float(1.613056249016E-14)
  ["VIII"]=>
  float(3.3395547440986E-28)
  ["IX"]=>
  float(9.4198561715257E-42)
  ["X"]=>
  float(2.5840062509191E-7)
  ["XI"]=>
  float(4.6985969966915E-21)
  ["XII"]=>
  float(1.6124316621054E-34)
  ["XIIA"]=>
  float(6.6577316325903E-48)
}
array(3) {
  ["degrees"]=>
  int(52)
  ["minutes"]=>
  int(39)
  ["seconds"]=>
  float(27.253085466713)
}
array(3) {
  ["degrees"]=>
  int(1)
  ["minutes"]=>
  int(43)
  ["seconds"]=>
  float(4.5176915241)
}