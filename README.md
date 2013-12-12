# geospatial - PHP Geospatial Extension

PHP Extension to handle common geospatial functions. The extension currently has implementations of the haversine and vincety folrmulas as well as a helmert transfomation function.

## Instalation
----------------------------

    git clone git@github.com:php-geospatial/geospatial.git
    cd geospatial
    phpize
    ./configure --enable-geospatial
    make
    sudo make install

    Then add the extension to an ini file e.g. /etc/php.ini

    extension = geospatial.so

## Usage

The extension makes use of the geojson standard format for co-ordinates. One important thing to note about this format is that points are specied longitude **first**.

e.g.
$greenwichObservatory = array(
	'type' => 'Point',
	'coordinates' => array( -0.001475 , 51.477811)
);


### Haversine

### Vincety

$flinders = array(
	'type' => 'Point',
	'coordinates' => array( $flindersPeakLong, $flindersPeakLat )
);
$buninyong = array(
	'type' => 'Point',
	'coordinates' => array( $buninyongLong, $buninyongLat )
);
var_dump(vincenty($flinders, $buninyong));


### Helmert Transformation

