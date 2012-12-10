--TEST--
Great Yarmouth Example From OS docs to OS Grid Eastings and Northings
--FILE--
<?php

var_dump(coord_to_eastings_northings(52.65757030556, 1.7179215833, true));
?>
--EXPECT--
array(15) {
  ["eastings"]=>
  float(651409.90293591)
  ["northings"]=>
  float(313177.27034344)
  ["e2"]=>
  float(0.0066705397615973)
  ["aF0"]=>
  float(6375020.480989)
  ["nu"]=>
  float(6388502.3332725)
  ["rho"]=>
  float(6372756.4398838)
  ["eta2"]=>
  float(0.0024708136168787)
  ["M"]=>
  float(406688.29596064)
  ["I"]=>
  float(306688.29596064)
  ["II"]=>
  float(1540407.9092342)
  ["III"]=>
  float(156068.75423782)
  ["IIIA"]=>
  float(-20671.123010934)
  ["IV"]=>
  float(3875120.5748557)
  ["V"]=>
  float(-170000.78207806)
  ["VI"]=>
  float(-101344.70431927)
}