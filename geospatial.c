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

zend_class_entry *php_geospatial_fc_entry;

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

ZEND_BEGIN_ARG_INFO_EX(coord_to_eastings_northings_args, 0, 0, 3)
	ZEND_ARG_INFO(0, latitude)
	ZEND_ARG_INFO(0, longitude)
	ZEND_ARG_INFO(0, returnIntermediateValues)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(os_grid_letters_args, 0, 0, 2)
	ZEND_ARG_INFO(0, eastings)
	ZEND_ARG_INFO(0, northings)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(os_grid_numeric_args, 0, 0, 2)
	ZEND_ARG_INFO(0, eastings)
	ZEND_ARG_INFO(0, northings)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(eastings_northings_to_coords_args, 0, 0, 3)
	ZEND_ARG_INFO(0, eastings)
	ZEND_ARG_INFO(0, northings)
	ZEND_ARG_INFO(0, returnIntermediateValues)
ZEND_END_ARG_INFO()

/* {{{ geospatial_functions[]
 *
 * Every user visible function must have an entry in geospatial_functions[].
 */
#define PHP_GEOSPATIAL_FC_NAME "LatLong"
const zend_function_entry geospatial_functions[] = {
	PHP_ME(LatLong, __construct, NULL, ZEND_ACC_PUBLIC|ZEND_ACC_CTOR)
	PHP_ME(LatLong, getHaversineDistance, NULL, ZEND_ACC_PUBLIC)
	PHP_ME(LatLong, transformDatum, NULL, ZEND_ACC_PUBLIC)
	PHP_FE(haversine, haversine_args)
	PHP_FE(helmert, helmert_args)
	PHP_FE(polar_to_cartesian, polar_to_cartesian_args)
	PHP_FE(cartesian_to_polar, cartesian_to_polar_args)
	PHP_FE(transform_datum, transform_datum_args)
	PHP_FE(dms_to_decimal, dms_to_decimal_args)
	PHP_FE(decimal_to_dms, decimal_to_dms_args)
	PHP_FE(coord_to_eastings_northings, coord_to_eastings_northings_args)
	PHP_FE(eastings_northings_to_coords, eastings_northings_to_coords_args)
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
	 zend_class_entry ce;
      INIT_CLASS_ENTRY(ce, PHP_GEOSPATIAL_FC_NAME,
                       geospatial_functions);
      php_geospatial_fc_entry =
            zend_register_internal_class(&ce TSRMLS_CC);
      // php_geospatial_fc_entry->create_object = create_ellipsoid;
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

zend_object_value create_ellipsoid(zend_class_entry *class_type TSRMLS_DC) {
  zend_object_value retval;
  geo_ellipsoid *intern;
  zval *tmp;

  // allocate the struct we're going to use
  intern = (geo_ellipsoid*)emalloc(sizeof(geo_ellipsoid));
  memset(intern, 0, sizeof(geo_ellipsoid));

  // // create a table for class properties
  zend_object_std_init(&intern->std, class_type TSRMLS_CC);

  // // create a destructor for this struct
  retval.handle = zend_objects_store_put(intern, (zend_objects_store_dtor_t) zend_objects_destroy_object, free_ellipsoid, NULL TSRMLS_CC);
  retval.handlers = zend_get_std_object_handlers();

  return retval;
}

void free_ellipsoid(void *object TSRMLS_DC) {

  geo_ellipsoid *ellipsoid = (geo_ellipsoid*)object;
  zend_hash_destroy(ellipsoid->std.properties);
  FREE_HASHTABLE(ellipsoid->std.properties);

  efree(ellipsoid);
}


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

geo_cartesian geospatial_helmert(double x, double y, double z, geo_helmert_constants helmert_consts)
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

double geospatial_meridional_arc(double Phi) {
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
	double degrees, minutes;
	double seconds, decimal;
	char *direction = "";
	int direction_len, sign;

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
	double decimal, seconds, absoluteDecimal;
	int degrees, minutes;
	char *direction;
	char *coordinate;
	int coordinate_len = 0;

	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "d|s", &decimal, &coordinate, &coordinate_len) == FAILURE) {
		return;
	}
	absoluteDecimal = fabs(decimal);
	if (coordinate_len > 0) {
		if (strcmp(coordinate, "longitude") == 0) {
			direction = decimal < 1 ? "W" : "E";
		} else {
			direction = decimal < 1 ? "S" : "N";
		}
		degrees = (int) absoluteDecimal;
	} else {
		degrees = (int) decimal;
	}

	array_init(return_value);
	minutes = absoluteDecimal * 60 - abs(degrees) * 60;
	seconds = absoluteDecimal * 3600 - abs(degrees) * 3600 - minutes * 60;
	add_assoc_long(return_value, "degrees", degrees);
	add_assoc_long(return_value, "minutes", minutes);
	add_assoc_double(return_value, "seconds", seconds);
	if (coordinate_len > 0) {
		add_assoc_string(return_value, "direction", direction, 1);
	}
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
	point = geospatial_helmert(x, y, z, helmert_constants);
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
		sprintf(gridLetters, "%c%c%04d%04d", letters[firstKey], letters[secondKey], remainderEast, remainderNorth);
		RETVAL_STRINGL(gridLetters, 10, 1);
	}	else {
		remainderEast = ((int) eastings ) / 100 % 1000;
		remainderNorth = ((int) northings ) /100 % 1000;
		sprintf(gridLetters, "%c%c%03d%03d", letters[firstKey], letters[secondKey], remainderEast, remainderNorth);
		RETVAL_STRINGL(gridLetters, 8, 1);
	}
	return;
}

PHP_FUNCTION(os_grid_numeric)
{
	double eastings, northings;
	char east[5], north[5];
	int length;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dd", &eastings, &northings) == FAILURE) {
		return;
	}
	array_init(return_value);
	sprintf(north, "%05d", (int) (northings / 100));
	sprintf(east, "%05d", (int) (eastings / 100));
	add_next_index_stringl(return_value, east, 5, 1);
	add_next_index_stringl(return_value, north, 5, 1);
}

geo_eta2 eta_sq(double Phi)
{
	double aF0, e2, nu, rho, eta2, sinPhi, sinPhi2;
	sinPhi = sin(Phi);
	e2 = (pow(airy_1830.a ,2) - pow(airy_1830.b, 2)) / pow(airy_1830.a, 2);
	sinPhi2 = pow(sinPhi, 2);
	aF0 = airy_1830.a * F_0;
	nu = aF0 / sqrt(1 - (e2 * sinPhi2));
	rho = (nu * (1 - e2))/(1 - e2 * sinPhi2);
	eta2 = nu / rho - 1;
	geo_eta2 etaSq = {e2, aF0, nu, rho, eta2};
	return etaSq;
}
coords_calculation eastings_northings_to_coords(double eastings, double northings)
{
	double latitude, longitude, phiDerivative, difference;
	double tanPhi, cosPhi, tanPhi2, aF0, M;
	coords_calculation p;
	geo_eta2 eta2;

	aF0 = airy_1830.a * F_0;
	difference = eastings - E_0;
	M = 0;
	phiDerivative = PHI_0;
	while ((northings - N_0 - M) >= 0.0001) {
		phiDerivative = (northings - N_0 - M) / aF0 + phiDerivative;
		M = geospatial_meridional_arc(phiDerivative);
	}
	cosPhi = cos(phiDerivative);
	tanPhi = tan(phiDerivative);
	tanPhi2 = pow(tanPhi, 2);
	eta2 = eta_sq(phiDerivative);

	p.phiDerivative = phiDerivative;
	p.eta2 = eta2;
	p.M = M;

	p.VII = tanPhi / (2 * eta2.rho * eta2.nu);
	p.VIII = tanPhi  / (24 * eta2.rho * pow(eta2.nu, 3)) * (5 + 3 * tanPhi2 + eta2.eta2 - 9 * tanPhi2 * eta2.eta2);
	p.IX = tanPhi / (720 * eta2.rho * pow(eta2.nu, 5)) * (61 + 90 * tanPhi2 + 45 * pow(tanPhi2, 2));
	p.X = 1 / (cosPhi * eta2.nu);
	p.XI = 1 / (cosPhi * 6 * pow(eta2.nu, 3)) * (eta2.nu / eta2.rho + 2 * tanPhi2);
	p.XII = 1 / (cosPhi * 120 * pow(eta2.nu, 5)) * (5 + 28 * tanPhi2 + 24 * pow(tanPhi2, 2));
	p.XIIA = 1 / (cosPhi * 5040 * pow(eta2.nu, 7)) * (61 + 662 * tanPhi2 + 1320 * pow(tanPhi2, 2) + 720 * pow(tanPhi ,6));
	p.latitude = phiDerivative - p.VII * pow(difference, 2) + p.VIII * pow(difference, 4) - p.IX * pow(difference, 6);
	p.longitude = LAMDBA_0 + p.X * difference - p.XI * pow(difference, 3) + p.XII * pow(difference, 5) - p.XIIA * pow(difference, 7);
	return p;
}


PHP_FUNCTION(eastings_northings_to_coords)
{
	double eastings, northings;
	geo_eta2 eta2;
	bool returnIntermediateValues;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dd|b", &eastings, &northings, &returnIntermediateValues) == FAILURE) {
		return;
	}

	coords_calculation coords = eastings_northings_to_coords(eastings, northings);
	eta2 = coords.eta2;

	array_init(return_value);
	add_assoc_double(return_value, "latitude", coords.latitude  / GEO_DEG_TO_RAD);
	add_assoc_double(return_value, "longitude", coords.longitude  / GEO_DEG_TO_RAD);
	if (returnIntermediateValues) {
		add_assoc_double(return_value, "phiDerivative", coords.phiDerivative);
		add_assoc_double(return_value, "M", coords.M);
		add_assoc_double(return_value, "nu", eta2.nu);
		add_assoc_double(return_value, "rho", eta2.rho);
		add_assoc_double(return_value, "eta2", eta2.eta2);
		add_assoc_double(return_value, "VII", coords.VII);
		add_assoc_double(return_value, "VIII", coords.VIII);
		add_assoc_double(return_value, "IX", coords.IX);
		add_assoc_double(return_value, "X", coords.X);
		add_assoc_double(return_value, "XI", coords.XI);
		add_assoc_double(return_value, "XII", coords.XII);
		add_assoc_double(return_value, "XIIA", coords.XIIA);
	}
}

eastings_northings_calculation coord_to_eastings_northings(double Phi, double Lambda)
{

	double n, bF0, sinPhi, cosPhi, tanPhi2;
	double rho, M, I, II, III, IIIA,IV, V, VI, P;
	geo_eta2 eta2;

	sinPhi = sin(Phi);
	cosPhi = cos(Phi);
	tanPhi2 = pow(tan(Phi), 2);

	eastings_northings_calculation results;
	eta2 = eta_sq(Phi);
	results.M = geospatial_meridional_arc(Phi);
	results.I = results.M + N_0;
	results.II = eta2.nu / 2 * sinPhi * cosPhi;
	results.III = eta2.nu / 24 * sinPhi * pow(cosPhi, 3) * (5 - tanPhi2 + 9 * eta2.eta2);
	results.IIIA = eta2.nu / 720 * sinPhi * pow(cosPhi, 5) * (61 - 58 * tanPhi2 + pow(tanPhi2, 2));
	results.IV = eta2.nu * cos(Phi);
	results.V = eta2.nu / 6 * pow(cosPhi, 3) * (eta2.nu / eta2.rho - tanPhi2);
	results.VI = eta2.nu / 120 * pow(cosPhi, 5) * (5 - 18 * tanPhi2 + pow(tanPhi2, 2) + 14 * eta2.eta2 -58 * tanPhi2 * eta2.eta2);
	results.eta2 = eta2;
	P = Lambda - LAMDBA_0;
	results.northings = results.I + results.II * pow(P, 2) + results.III * pow(P, 4) + results.IIIA * pow(P, 6);
	results.eastings =  E_0 + results.IV * P + results.V * pow(P, 3) + results.VI * pow(P, 5);
	return results;
}

PHP_FUNCTION(coord_to_eastings_northings)
{
	double latitude, longitude;
	double Phi, Lambda;
	bool returnIntermediateValues;
	geo_eta2 eta2;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dd|b", &latitude, &longitude, &returnIntermediateValues) == FAILURE) {
		return;
	}
	Phi = latitude * GEO_DEG_TO_RAD;
	Lambda = longitude * GEO_DEG_TO_RAD;
	eastings_northings_calculation east_north = coord_to_eastings_northings(Phi, Lambda);
	eta2 = east_north.eta2;
	array_init(return_value);

	add_assoc_double(return_value, "eastings", east_north.eastings);
	add_assoc_double(return_value, "northings", east_north.northings);
	if (returnIntermediateValues) {
		add_assoc_double(return_value, "e2", eta2.e2);
		add_assoc_double(return_value, "aF0", eta2.aF0);
		add_assoc_double(return_value, "nu", eta2.nu);
		add_assoc_double(return_value, "rho", eta2.rho);
		add_assoc_double(return_value, "eta2", eta2.eta2);
		add_assoc_double(return_value, "M", east_north.M);
		add_assoc_double(return_value, "I", east_north.I);
		add_assoc_double(return_value, "II", east_north.II);
		add_assoc_double(return_value, "III", east_north.III);
		add_assoc_double(return_value, "IIIA", east_north.IIIA);
		add_assoc_double(return_value, "IV", east_north.IV);
		add_assoc_double(return_value, "V", east_north.V);
		add_assoc_double(return_value, "VI", east_north.VI);
	}
}

geo_lat_long geospatial_transform_datum(double latitude, double longitude,
						 int from_reference_ellipsoid, int to_reference_ellipsoid)
{
	geo_cartesian point, converted_point;

	geo_ellipsoid eli_from = get_ellipsoid(from_reference_ellipsoid);
	geo_ellipsoid eli_to = get_ellipsoid(to_reference_ellipsoid);
	point = polar_to_cartesian(latitude, longitude, eli_from);
	geo_helmert_constants helmert_constants = get_helmert_constants(from_reference_ellipsoid, to_reference_ellipsoid);
	converted_point = geospatial_helmert(point.x, point.y, point.z, helmert_constants);
	return cartesian_to_polar(converted_point.x, converted_point.y, converted_point.z, eli_to);
}

/* {{{ proto transform_datum(double latitude, double longitude, long from_reference_ellipsoid, long to_reference_ellipsoid)
 * Unified function to transform projection of geo-cordinates between datums */
PHP_FUNCTION(transform_datum)
{
	double latitude, longitude;
	long from_reference_ellipsoid, to_reference_ellipsoid;
	geo_lat_long polar;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ddll", &latitude, &longitude, &from_reference_ellipsoid, &to_reference_ellipsoid) == FAILURE) {
		return;
	}

	polar = geospatial_transform_datum(latitude, longitude, from_reference_ellipsoid, to_reference_ellipsoid);

	array_init(return_value);
	add_assoc_double(return_value, "lat", polar.latitude);
	add_assoc_double(return_value, "long", polar.longitude);
	add_assoc_double(return_value, "height", polar.height);
}
/* }}} */




int geospatial_get_sign_from_direction(char *direction)
{
	return strcmp(direction, "S") == 0 || strcmp(direction, "W") == 0 ? -1 : 1;
}

PHP_METHOD(LatLong, __construct)
{
    char *lat_str = NULL;
    int lat_str_len = 0;
    char *long_str = NULL;
    int long_str_len = 0, sign;
    double latitude, longitude;
    char lat_dir[2], long_dir[2];
    double minutes;
    int degrees;
    long reference_ellipsoid = GEO_WGS84;
	zval *objvar = getThis();

	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ss|l", &lat_str, &lat_str_len, &long_str, &long_str_len, &reference_ellipsoid) == FAILURE) {
		return;
	}

	if (!objvar) {
		php_error_docref(NULL TSRMLS_CC, E_WARNING, "Constructor called statically!");
		RETURN_FALSE;
	}


	if (sscanf(lat_str, "%lf %c", &latitude, lat_dir)) {
		sign = geospatial_get_sign_from_direction(lat_dir);
	    add_property_double(objvar, "latitude", sign*latitude);
	} else if (sscanf(lat_str, "%d° %lf\" %c", &degrees, &minutes, lat_dir)) {
		sign = geospatial_get_sign_from_direction(lat_dir);
	    add_property_double(objvar, "latitude", sign*degrees + minutes / 60.0);
	}
	if (sscanf(long_str, "%lf %c", &longitude, long_dir)) {
		sign = geospatial_get_sign_from_direction(long_dir);
	    add_property_double(objvar, "longitude", sign*longitude);
	}

	add_property_long(objvar, "reference_ellipsoid", reference_ellipsoid);
}

double geospatial_haversine(double fromLat, double fromLong, double toLat, double toLong, double radius)
{
	double deltaLat, deltaLong,  latH, longH, result;

	deltaLat = (fromLat - toLat) * GEO_DEG_TO_RAD;
	deltaLong = (fromLong - toLong) * GEO_DEG_TO_RAD;

	latH = pow(sin(deltaLat * 0.5), 2);
	longH = pow(sin(deltaLong * 0.5), 2);

	result = cos(fromLat * GEO_DEG_TO_RAD) * cos(toLat * GEO_DEG_TO_RAD);
	return radius * 2.0 * asin(sqrt(latH + result * longH));
}

/* {{{ proto haversine(double fromLat, double fromLong, double toLat, double toLong [, double radius ])
 * Calculates the greater circle distance between the two lattitude/longitude pairs */
PHP_FUNCTION(haversine)
{
	double fromLat, fromLong, toLat, toLong;
	double radius = GEO_EARTH_RADIUS, result;
	if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dddd|d", &fromLat, &fromLong, &toLat, &toLong, &radius) == FAILURE) {
		return;
	}

	result = geospatial_haversine(fromLat, fromLong, toLat, toLong, radius);
	RETURN_DOUBLE(result);
}
/* }}} */

PHP_METHOD(LatLong, transformDatum)
{
	long  to_reference_ellipsoid = GEO_WGS84;
    zval *thisobjvar = getThis();
	geo_cartesian point, converted_point;
	geo_lat_long polar;
    if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "l",
								          &to_reference_ellipsoid) == FAILURE) {
        RETURN_NULL();
	}

	zval **fromLat;
	if (zend_hash_find(Z_OBJPROP_P(thisobjvar),
					           "latitude", sizeof("latitude"), (void **)&fromLat) == FAILURE) {
		    return;
	}
	zval **fromLong;
	if (zend_hash_find(Z_OBJPROP_P(thisobjvar),
	                           "longitude", sizeof("longitude"), (void**)&fromLong) == FAILURE) {
	         return;
	}
	zval **fromEllipsoid;
	if (zend_hash_find(Z_OBJPROP_P(thisobjvar),
	                           "reference_ellipsoid", sizeof("reference_ellipsoid"), (void**)&fromEllipsoid) == FAILURE) {
	         return;
	}


	polar = geospatial_transform_datum(Z_DVAL_PP(fromLat), Z_DVAL_PP(fromLong), Z_LVAL_PP(fromEllipsoid), to_reference_ellipsoid);

	array_init(return_value);
	add_assoc_double(return_value, "lat", polar.latitude);
	add_assoc_double(return_value, "long", polar.longitude);
	add_assoc_double(return_value, "height", polar.height);
}


PHP_METHOD(LatLong, getHaversineDistance)
{
	double deltaLat, deltaLong;
	double radius = GEO_EARTH_RADIUS, result;
	zval *objvar;
    zval *thisobjvar = getThis();
    if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "O|d",
								          &objvar, php_geospatial_fc_entry, &radius) == FAILURE) {
			        RETURN_NULL();
	}
	zval **fromLat;
	if (zend_hash_find(Z_OBJPROP_P(thisobjvar),
					           "latitude", sizeof("latitude"), (void **)&fromLat) == FAILURE) {
		    return;
	}
	zval **fromLong;
	if (zend_hash_find(Z_OBJPROP_P(thisobjvar),
	                           "longitude", sizeof("longitude"), (void**)&fromLong) == FAILURE) {
	         return;
	}

	zval **toLat;
	if (zend_hash_find(Z_OBJPROP_P(objvar),
	                            "latitude", sizeof("latitude"), (void**)&toLat) == FAILURE) {
	          return;
	}
	 zval **toLong;
	     if (zend_hash_find(Z_OBJPROP_P(objvar),
	                           "longitude", sizeof("longitude"), (void**)&toLong) == FAILURE) {
	           return;
	}

	result = geospatial_haversine(Z_DVAL_PP(fromLat), Z_DVAL_PP(fromLong), Z_DVAL_PP(toLat), Z_DVAL_PP(toLong), radius);
	RETURN_DOUBLE(result);
}

/*
 * Local variables:
 * tab-width: 4
 * c-basic-offset: 4
 * End:
 * vim600: noet sw=4 ts=4 fdm=marker
 * vim<600: noet sw=4 ts=4
 */
