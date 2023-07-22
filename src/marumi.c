/*
 * marumi.c Earth curvature by oue
 */

#include "marumi.h"
#include <math.h>


#define EARTH_RADIUS_km    6.37E3
#define EFFECTIVE_RADIUS_km  ((EARTH_RADIUS_km)*(1.33333))

void marumi( const double az_rad, const double el_rad, const double r_km,
	     struct MARUMI *mp )
{
  
  double c;
  
  c = r_km * cos( el_rad ) / ( EFFECTIVE_RADIUS_km + r_km * sin( el_rad ) );

  if(  fabs(c) < 1.0e-6 ) {
    mp->dx_km = 0.0;
    mp->dy_km = 0.0;
    mp->dz_km = r_km;
  }
  else{
    c = atan( c );
    c = el_rad + c;
    mp->dx_km = r_km * cos( c ) * sin( az_rad );
    mp->dy_km = r_km * cos( c ) * cos( az_rad );
    mp->dz_km = r_km * sin( c ) 
      - pow( r_km * cos(c), 2.0 ) / 2.0 / EFFECTIVE_RADIUS_km;
  }
  
  return;

}
