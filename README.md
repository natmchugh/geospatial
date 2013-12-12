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

The extension makes use of the geojson standard format for co-ordinates. One important thing to note about this format is that points are specied longitude, latitude.

### Haversine

### Vincety

### Helmert Transformation

