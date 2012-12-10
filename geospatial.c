/*
  +----------------------------------------------------------------------+
  | PHP Version 5                                                        |
  +----------------------------------------------------------------------+
  | Copyright (c) 1997-2012 The PHP Group                                |
  +----------------------------------------------------------------------+
  | This source file is subject to version 3.01 of the PHP license,      |
  | that is bundled with this package in the file LICENSE, and is        |
  | available through the world-wide-web at the following url:           |
  | http://www.php.net/license/3_01.txt                                  |
  | If you did not receive a copy of the PHP license and are unable to   |
  | obtain it through the world-wide-web, please send a note to          |
  | license@php.net so we can mail you a copy immediately.               |
  +----------------------------------------------------------------------+
  | Author:                                                              |
  +----------------------------------------------------------------------+
*/


/* $Id$ */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "php.h"
#include "php_ini.h"
#include "ext/standard/info.h"
#include "php_geospatial.h"

ZEND_BEGIN_ARG_INFO_EX(haversine_args, 0, 0, 4)
	ZEND_ARG_INFO(0, fromLatitude)
	ZEND_ARG_INFO(0, fromLongitude)
	ZEND_ARG_INFO(0, toLatitude)
	ZEND_ARG_INFO(0, toLongitude)
	ZEND_ARG_INFO(0, radius)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(helmert_args, 0, 0, 3)
	ZEND_ARG_INFO(0, x)
	ZEND_ARG_INFO(0, y)
	ZEND_ARG_INFO(0, z)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(polar_to_cartesian_args, 0, 0, 3)
	ZEND_ARG_INFO(0, latitude)
	ZEND_ARG_INFO(0, longitude)
	ZEND_ARG_INFO(0, reference_ellipsoid)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(cartesian_to_polar_args, 0, 0, 4)
	ZEND_ARG_INFO(0, x)
	ZEND_ARG_INFO(0, y)
	ZEND_ARG_INFO(0, z)
	ZEND_ARG_INFO(0, reference_ellipsoid)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(transform_datum_args, 0, 0, 4)
	ZEND_ARG_INFO(0, latitude)
	ZEND_ARG_INFO(0, longitude)
	ZEND_ARG_INFO(0, from_reference_ellipsoid)
	ZEND_ARG_INFO(0, to_reference_ellipsoid)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(dms_to_decimal_args, 0, 0, 4)
	ZEND_ARG_INFO(0, degrees)
	ZEND_ARG_INFO(0, minutes)
	ZEND_ARG_INFO(0, seconds)
	ZEND_ARG_INFO(0, direction)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(decimal_to_dms_args, 0, 0, 2)
	ZEND_ARG_INFO(0, decimal)
	ZEND_ARG_INFO(0, coordinate)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(coord_to_eastings_northings_args, 0, 0, 2)
	ZEND_ARG_INFO(0, latitude)
	ZEND_ARG_INFO(0, longitude)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(os_grid_letters_args, 0, 0, 2)
	ZEND_ARG_INFO(0, eastings)
	ZEND_ARG_INFO(0, northings)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(os_grid_numeric_args, 0, 0, 2)
	ZEND_ARG_INFO(0, eastings)
	ZEND_ARG_INFO(0, northings)
ZEND_END_ARG_INFO()


/* {{{ geospatial_functions[]
 *
 * Every user visible function must have an entry in geospatial_functions[].
 */
const zend_function_entry geospatial_functions[] = {
	PHP_FE(haversine, haversine_args)
	PHP_FE(helmert, helmert_args)
	PHP_FE(polar_to_cartesian, polar_to_cartesian_args)
	PHP_FE(cartesian_to_polar, cartesian_to_polar_args)
	PHP_FE(transform_datum, transform_datum_args)
	PHP_FE(dms_to_decimal, dms_to_decimal_args)
	PHP_FE(decimal_to_dms, decimal_to_dms_args)
	PHP_FE(coord_to_eastings_northings, coord_to_eastings_northings_args)
	PHP_FE(os_grid_letters, os_grid_letters_args)
	PHP_FE(os_grid_numeric, os_grid_numeric_args)
	/* End of functions */
	{ NULL, NULL, NULL }
};
/* }}} */

/* {{{ geospatial_module_entry
 */
zend_module_entry geospatial_module_entry = {
#if ZEND_MODULE_API_NO >= 20010901
	STANDARD_MODULE_HEADER,
#endif
	"geospatial",
	geospatial_functions,
	PHP_MINIT(geospatial),
	NULL,
	NULL,
	NULL,
	PHP_MINFO(geospatial),
#if ZEND_MODULE_API_NO >= 20010901
	"0.1", /* Replace with version number for your extension */
#endif
	STANDARD_MODULE_PROPERTIES
};
/* }}} */

#ifdef COMPILE_DL_GEOSPATIAL
ZEND_GET_MODULE(geospatial)
#endif

/* {{{ PHP_MINIT_FUNCTION
 */
PHP_MINIT_FUNCTION(geospatial)
{
	REGISTER_DOUBLE_CONSTANT("GEO_DEG_TO_RAD", GEO_DEG_TO_RAD, CONST_CS | CONST_PERSISTENT);
	REGISTER_DOUBLE_CONSTANT("GEO_EARTH_RADIUS", GEO_EARTH_RADIUS, CONST_CS | CONST_PERSISTENT);
	REGISTER_LONG_CONSTANT("GEO_AIRY_1830", GEO_AIRY_1830, CONST_CS | CONST_PERSISTENT);
	REGISTER_LONG_CONSTANT("GEO_WGS84", GEO_WGS84, CONST_CS | CONST_PERSISTENT);
	return SUCCESS;
}
/* }}} */

/* {{{ PHP_MSHUTDOWN_FUNCTION
 */
PHP_MSHUTDOWN_FUNCTION(geospatial)
{
	return SUCCESS;
}
/* }}} */

/* {{{ PHP_MINFO_FUNCTION
 */
PHP_MINFO_FUNCTION(geospatial)
{
	php_info_print_table_start();
	php_info_print_table_header(2, "Geospatial functions", "enabled");
	php_info_print_table_end();
}
/* }}} */

geo_ellipsoid get_ellipsoid(long ellipsoid_const)
{
	switch (ellipsoid_const) {
		case GEO_AIRY_1830:
			return airy_1830;
			break;
		case GEO_WGS84:
		default:
			return wgs84;
			break;
	}
}

geo_helmert_constants get_helmert_constants(long from, long to)
{
	switch (from - to) {
		case 1:
			return osgb36_wgs84;
			break;
		default:
		case -1:
			return wgs84_osgb36;
			break;
	}
}

geo_cartesian helmert(double x, double y, double z, geo_helmert_constants helmert_consts)
{
	double rX, rY, rZ;
	double xOut, yOut, zOut;
	double scale_change;
	geo_cartesian point;
	scale_change = 1 + (helmert_consts.scale_change);
	rX = helmert_consts.rotation_x / GEO_SEC_IN_DEG * GEO_DEG_TO_RAD;
	rY = helmert_consts.rotation_y / GEO_SEC_IN_DEG * GEO_DEG_TO_RAD;
	rZ = helmert_consts.rotation_z / GEO_SEC_IN_DEG * GEO_DEG_TO_RAD;

	xOut = helmert_consts.translation_x + ((x - (rZ * y) + (rY * z)) * scale_change);

	yOut =  helmert_consts.translation_y + (((rZ * x) + y - (rX * z)) * scale_change);

	zOut = helmert_consts.translation_z + (((-1 * rY * x) + (rX * y) + z) * scale_change);

	point.x = xOut;
	point.y = yOut;
	point.z = zOut;
	return point;
}

geo_cartesian polar_to_cartesian(double latitude, double longitude, geo_ellipsoid eli)
{
	double x, y, z;

	geo_cartesian point;
	double phi = latitude * GEO_DEG_TO_RAD;
	double lambda = longitude * GEO_DEG_TO_RAD;
	double eSq = ((eli.a * eli.a)  - (eli.b * eli.b))  /  (eli.a * eli.a);
	double nu = eli.a / sqrt(1 - (eSq * sin(phi) * sin(phi)));
	x = nu + HEIGHT;
	x *= cos(phi) * cos(lambda);
	y = nu + HEIGHT;
	y *= cos(phi) * sin(lambda);
	z = ((1 - eSq) * nu) + HEIGHT;
	z*= sin(phi);
	point.x = x;
	point.y = y;
	point.z = z;
	return point;
}

double meridional_arc(double Phi) {
	double sum, difference;
	double aF0, bF0, n;
	double A, B, C, D;

	aF0 = airy_1830.a * F_0;
	bF0 = airy_1830.b * F_0;
	n = (aF0 - bF0) / (aF0 + bF0);
	sum = Phi + PHI_0;
	difference = Phi - PHI_0;

	A = (1 + n + 5.0 / 4.0 * pow(n, 2) + 5.0 / 4.0 * pow(n, 3)) * difference;
	B = (3 * n + 3 * pow(n, 2)  + 21.0 / 8.0 * pow(n, 3))  * sin(difference) * cos(sum);
	C = (15.0 / 8.0 * pow(n, 2) + 15.0 / 8.0 * pow(n, 3)) * sin(2 * difference) * cos(2 * sum);
	D = 35.0 / 24.0 * pow(n, 3) * sin(3 * difference) * cos(3 * sum);
	return bF0 * (A - B + C - D);
}


geo_lat_long cartesian_to_polar(double x, double y, double z, geo_ellipsoid eli)
{

	double latitude, longitude;
	double nu, lambda, h;
	geo_lat_long polar;

	/* aiming for 1m accuracy */
	double precision = 0.1 / eli.a;
	double eSq = ((eli.a * eli.a)  - (eli.b * eli.b))  /  (eli.a * eli.a);
	double p = sqrt(x * x + y * y);
	double phi = atan2(z, p * (1 - eSq));
	double phiP = 2 * M_PI;

	while (abs(phi - phiP) > precision) {
		nu = eli.a / sqrt(1 - eSq * sin(phi) * sin(phi));
		phiP = phi;
		phi = atan2(z + eSq * nu * sin(phi), p);
	}

	lambda = atan2(y ,x);
	h = p / cos(phi) - nu;
	polar.latitude = phi / GEO_DEG_TO_RAD;
	polar.longitude = lambda / GEO_DEG_TO_RAD;
	polar.height = h;

	return polar;
}

/* {{{ proto dms_to_decimal(double degrees, double minutes, double seconds [,string direction])
 * Convert degrees, minutes & seconds values to decimal degrees */
PHP_FUNCTION(dms_to_decimal)
{
	double degrees, minutes, sign;
	double seconds, decimal;
	char *direction = "";
	int direction_len;

	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ddd|s", &degrees, &minutes, &seconds, &direction, &direction_len) == FAILURE) {
		return;
	}

	if (strcmp("", direction) == 0) {
		sign = degrees > 1 ? 1 : -1;
	} else {
		sign = strcmp(direction, "S") == 0 || strcmp(direction, "W") == 0 ? -1 : 1;
	}

	decimal = abs(degrees) + minutes / 60 + seconds / 3600;
	decimal *= sign;
	RETURN_DOUBLE(decimal);
}
/* }}} */

/* {{{ proto decimal_to_dms(double decimal, string coordinate)
 * Convert decimal degrees value to whole degrees and minutes and decimal seconds */
PHP_FUNCTION(decimal_to_dms)
{
	double decimal, seconds;
	int degrees, minutes;
	char *direction;
	char *coordinate;
	int coordinate_len;

	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ds", &decimal, &coordinate, &coordinate_len) == FAILURE) {
		return;
	}

	if (strcmp(coordinate, "longitude") == 0) {
		direction = decimal < 1 ? "W" : "E";
	} else {
		direction = decimal < 1 ? "S" : "N";
	}

	array_init(return_value);
	decimal = fabs(decimal);
	degrees = (int) decimal;
	minutes = decimal * 60 - degrees * 60;
	seconds = decimal * 3600 - degrees * 3600 - minutes * 60;
	add_assoc_long(return_value, "degrees", degrees);
	add_assoc_long(return_value, "minutes", minutes);
	add_assoc_double(return_value, "seconds", seconds);
	add_assoc_string(return_value, "direction", direction, 1);
}
/* }}} */

/* {{{ proto helmert(double x, double y, double z [, long from_reference_ellipsoid, long to_reference_ellipsoid])
 * Convert polar ones (latitude, longitude) tp cartesian co-ordiantes (x, y, z)  */
PHP_FUNCTION(helmert)
{
	double x, y, z;
	geo_cartesian point;
	long from_reference_ellipsoid, to_reference_ellipsoid;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ddd|ll", &x, &y, &z, &from_reference_ellipsoid, &to_reference_ellipsoid) == FAILURE) {
		return;
	}

	array_init(return_value);
	geo_helmert_constants helmert_constants = get_helmert_constants(from_reference_ellipsoid, to_reference_ellipsoid);
	point = helmert(x, y, z, helmert_constants);
	add_assoc_double(return_value, "x", point.x);
	add_assoc_double(return_value, "y", point.y);
	add_assoc_double(return_value, "z", point.z);
}
/* }}} */

/* {{{ proto polar_to_cartesian(double latitude, double longitude[, long reference_ellipsoid])
 * Convert polar ones (latitude, longitude) tp cartesian co-ordiantes (x, y, z)  */
PHP_FUNCTION(polar_to_cartesian)
{
	double latitude, longitude;
	long reference_ellipsoid;
	geo_cartesian point;

	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dd|l", &latitude, &longitude, &reference_ellipsoid) == FAILURE) {
		return;
	}

	geo_ellipsoid eli = get_ellipsoid(reference_ellipsoid);
	array_init(return_value);
	point = polar_to_cartesian(latitude, longitude, eli);
	add_assoc_double(return_value, "x", point.x);
	add_assoc_double(return_value, "y", point.y);
	add_assoc_double(return_value, "z", point.z);
}
/* }}} */

/* {{{ proto cartesian_to_polar(double x, double y, double z [, long reference_ellipsoid])
 * Convert cartesian co-ordiantes (x, y, z) to polar ones (latitude, longitude) */
PHP_FUNCTION(cartesian_to_polar)
{
	double x, y, z;
	long reference_ellipsoid;
	geo_lat_long polar;

	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ddd|l", &x, &y, &z, &reference_ellipsoid) == FAILURE) {
		return;
	}

	geo_ellipsoid eli = get_ellipsoid(reference_ellipsoid);
	array_init(return_value);
	polar = cartesian_to_polar(x, y, z, eli);
	add_assoc_double(return_value, "lat", polar.latitude);
	add_assoc_double(return_value, "long", polar.longitude);
	add_assoc_double(return_value, "height", polar.height);
}
/* }}} */

PHP_FUNCTION(os_grid_letters)
{
	double eastings, northings;
	int hundredsKmNorth, hundredsKmEast, vertical500Count, horizontal500Count;
	int firstKey, secondKey, top, left, verticalCount, horizotalCount;
	int remainderEast, remainderNorth, numberFigures = 6;
	char letters[] = "ABCDEFGHJKLMNOPQRSTUVWXZ", gridLetters[2];
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dd|l", &eastings, &northings, &numberFigures) == FAILURE) {
		return;
	}
	hundredsKmNorth = floor(northings / 100000);
	vertical500Count = floor(hundredsKmNorth / 5);
	hundredsKmEast = floor(eastings / 100000);
	horizontal500Count = floor(hundredsKmEast / 5);
	firstKey = horizontal500Count +2 + 5* (3 - vertical500Count);
	top = 5*(vertical500Count+1);
	left = 5 * horizontal500Count;
	verticalCount = top - hundredsKmNorth - 1;
	horizotalCount = hundredsKmEast - left;
	secondKey = verticalCount * 5 + horizotalCount;
	if (8 == numberFigures) {
		remainderEast = ((int) eastings ) / 10 % 10000;
		remainderNorth = ((int) northings ) /10 % 10000;
		sprintf(gridLetters, "%c%c%4d%4d", letters[firstKey], letters[secondKey], remainderEast, remainderNorth);
		RETVAL_STRINGL(gridLetters, 10, 1);
	}	else {
		remainderEast = ((int) eastings ) / 100 % 1000;
		remainderNorth = ((int) northings ) /100 % 1000;
		sprintf(gridLetters, "%c%c%3d%3d", letters[firstKey], letters[secondKey], remainderEast, remainderNorth);
		RETVAL_STRINGL(gridLetters, 8, 1);
	}
	return;
}

PHP_FUNCTION(os_grid_numeric)
{
	double eastings, northings;
	char east[4], north[4];
	int length;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dd", &eastings, &northings) == FAILURE) {
		return;
	}
	array_init(return_value);
	sprintf(north, "%05d", ((int) northings / 100));
	sprintf(east, "%05d", ((int) eastings / 100));
	add_next_index_stringl(return_value, east, 5, 1);
	add_next_index_stringl(return_value, north, 5, 1);
}

PHP_FUNCTION(coord_to_eastings_northings)
{
	double latitude, longitude;
	double e2, n, nu, aF0, bF0, Phi, Lambda, sinPhi, cosPhi, sinPhi2, tanPhi2, eta2;
	double rho, M, I, II, III, IIIA,IV, V, VI, P;
	bool returnIntermediateValues;
	double eastings, northings;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dd|b", &latitude, &longitude, &returnIntermediateValues) == FAILURE) {
		return;
	}
	e2 = (pow(airy_1830.a ,2) - pow(airy_1830.b, 2)) / pow(airy_1830.a, 2);
	aF0 = airy_1830.a * F_0;
	Phi = latitude * GEO_DEG_TO_RAD;
	Lambda = longitude * GEO_DEG_TO_RAD;
	sinPhi = sin(Phi);
	cosPhi = cos(Phi);
	sinPhi2 = pow(sinPhi, 2);
	tanPhi2 = pow(tan(Phi), 2);
	nu = aF0 / sqrt(1 - (e2 * sinPhi2));
	rho = (nu * (1 - e2))/(1 - e2 * sinPhi2);
	eta2 = nu / rho - 1;
	M = meridional_arc(Phi);
	I = M + N_0;
	II = nu / 2 * sinPhi * cosPhi;
	III = nu / 24 * sinPhi * pow(cosPhi, 3) * (5 - tanPhi2 + 9 * eta2);
	IIIA = nu / 720 * sinPhi * pow(cosPhi, 5) * (61 - 58 * tanPhi2 + pow(tanPhi2, 2));
	IV = nu * cos(Phi);
	V = nu / 6 * pow(cosPhi, 3) * (nu / rho - tanPhi2);
	VI = nu / 120 * pow(cosPhi, 5) * (5 - 18 * tanPhi2 + pow(tanPhi2, 2) + 14 * eta2 -58 * tanPhi2 * eta2);
	P = Lambda - LAMDBA_0;
	northings = I + II * pow(P, 2) + III * pow(P, 4) + IIIA * pow(P, 6);
	eastings =  E_0 + IV * P + V * pow(P, 3) + VI * pow(P, 5);

	array_init(return_value);

	add_assoc_double(return_value, "eastings", eastings);
	add_assoc_double(return_value, "northings", northings);
	if (returnIntermediateValues) {
		add_assoc_double(return_value, "e2", e2);
		add_assoc_double(return_value, "aF0", aF0);
		add_assoc_double(return_value, "nu", nu);
		add_assoc_double(return_value, "rho", rho);
		add_assoc_double(return_value, "eta2", eta2);
		add_assoc_double(return_value, "M", M);
		add_assoc_double(return_value, "I", I);
		add_assoc_double(return_value, "II", II);
		add_assoc_double(return_value, "III", III);
		add_assoc_double(return_value, "IIIA", IIIA);
		add_assoc_double(return_value, "IV", IV);
		add_assoc_double(return_value, "V", V);
		add_assoc_double(return_value, "VI", VI);
	}
}


/* {{{ proto transform_datum(double latitude, double longitude, long from_reference_ellipsoid, long to_reference_ellipsoid)
 * Unified function to transform projection of geo-cordinates between datums */
PHP_FUNCTION(transform_datum)
{
	double latitude, longitude;
	long from_reference_ellipsoid, to_reference_ellipsoid;
	geo_cartesian point, converted_point;
	geo_lat_long polar;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ddll", &latitude, &longitude, &from_reference_ellipsoid, &to_reference_ellipsoid) == FAILURE) {
		return;
	}

	geo_ellipsoid eli_from = get_ellipsoid(from_reference_ellipsoid);
	geo_ellipsoid eli_to = get_ellipsoid(to_reference_ellipsoid);
	point = polar_to_cartesian(latitude, longitude, eli_from);
	geo_helmert_constants helmert_constants = get_helmert_constants(from_reference_ellipsoid, to_reference_ellipsoid);
	converted_point = helmert(point.x, point.y, point.z, helmert_constants);
	polar = cartesian_to_polar(converted_point.x, converted_point.y, converted_point.z, eli_to);

	array_init(return_value);
	add_assoc_double(return_value, "lat", polar.latitude);
	add_assoc_double(return_value, "long", polar.longitude);
	add_assoc_double(return_value, "height", polar.height);
}
/* }}} */

/* {{{ proto haversine(double fromLat, double fromLong, double toLat, double toLong [, double radius ])
 * Calculates the greater circle distance between the two lattitude/longitude pairs */
PHP_FUNCTION(haversine)
{
	double fromLat, fromLong, toLat, toLong, deltaLat, deltaLong;
	double radius = GEO_EARTH_RADIUS, latH, longH, result;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dddd|d", &fromLat, &fromLong, &toLat, &toLong, &radius) == FAILURE) {
		return;
	}

	deltaLat = (fromLat - toLat) * GEO_DEG_TO_RAD;
	deltaLong = (fromLong - toLong) * GEO_DEG_TO_RAD;

	latH = sin(deltaLat * 0.5);
	latH *= latH;
	longH = sin(deltaLong * 0.5);
	longH *= longH;

	result = cos(fromLat * GEO_DEG_TO_RAD) * cos(toLat * GEO_DEG_TO_RAD);
	result = radius * 2.0 * asin(sqrt(latH + result * longH));
	RETURN_DOUBLE(result);
}
/* }}} */

/*
 * Local variables:
 * tab-width: 4
 * c-basic-offset: 4
 * End:
 * vim600: noet sw=4 ts=4 fdm=marker
 * vim<600: noet sw=4 ts=4
 */
