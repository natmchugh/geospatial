# geospatial - PHP Geospatial Extension

PHP Extension to handle common geospatial functions. The extension currently has implementations of the haversine and vincenty's formulas as well as a helmert transfomation function.

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

The extension makes use of the geojson standard format for specying points in terms ofco-ordinates. One important thing to note about this format is that points are specied longitude **first** i.e. longitude, latitude.

e.g.

	$greenwichObservatory = array(
		'type' => 'Point',
		'coordinates' => array( -0.001475 , 51.477811)
	);


### Haversine

	$from = array(
		'type' => 'Point',
		'coordinates' => array( -104.88544, 39.06546 )
	);
	$to = array(
		'type' => 'Point',
		'coordinates' => array( -104.80, 39.06546 )
	);
	var_dump(haversine($to, $from));
	

### Vincenty's Formula

Vincenty's formula attempts to provide a more acurate distance between two points than the Haversine formula. Whereas the Haversine formula assumes a spherical earth the Vincenty method models the earth as an ellipsoid.

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

The Helmert transformation allows for the transfomation of points between different datums. It can for instance be used to convert between the WGS84 ellipsoid used by GPS systems and OSGB36 used by ordnance survey in the UK.

	$from = array('type' => 'Point', 'coordinates' => array( $long, $lat ) );
	
	$polar = transform_datum($from, GEO_WGS84, GEO_AIRY_1830);
	
	var_dump($polar);



