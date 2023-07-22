/*
 * marumi.h 
 */

#ifndef MSAKA_MARUMI_H_
#define MSAKA_MARUMI_H_


struct MARUMI {
  double dx_km, dy_km, dz_km;
}; 


void marumi( const double az_rad, const double el_rad, const double r_km,
	     struct MARUMI *mp );


#endif
