/*
freesteam - IAPWS-IF97 steam tables library
Copyright (C) 2004-2009  John Pye

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*//** @file
	Implementation details for backwards equations. Intended that
	this header file would only be used for internal tests etc.
*/
#ifndef FREESTEAM_BACKWARDS_IMPL_H
#define FREESTEAM_BACKWARDS_IMPL_H

#include "common.h"

#define SQ(X) ((X)*(X))
#define CUBE(X) ((X)*(X)*(X))

/* boundary between subregions 3a and 3b in region 3 for (p,h) */

#define REGION3_B3AB_PSTAR (1.e6)
#define REGION3_B3AB_HSTAR (1.e3)
#define REGION3_B3AB_ETA(H) ((H)/REGION3_B3AB_HSTAR)
#define REGION3_B3AB_PI(P) ((P)/REGION3_B3AB_PSTAR)

#define REGION3_B3AB_PH(P,H) (REGION3_B3AB_ETA(H) - (\
	2014.64004206875 \
	+ 3.74696550136983*REGION3_B3AB_PI(P) \
	- 2.19921901054187E-02 * SQ(REGION3_B3AB_PI(P)) \
	+ 8.7513168600995E-05 * CUBE(REGION3_B3AB_PI(P)) \
	))


/* boundary between subregions 2 and 3 in region 2 for (p,h) */

#define REGION2_B2BC_PSTAR (1.e6)
#define REGION2_B2BC_HSTAR (1.e3)
#define REGION2_B2BC_ETA(H) ((H)/REGION2_B2BC_HSTAR)
#define REGION2_B2BC_PI(P) ((P)/REGION2_B2BC_PSTAR)

#define REGION2_B2BC_PH(P,H) (\
	(REGION2_B2BC_PI(P) - (\
	905.84278514723 \
	- 0.67955786399241*REGION2_B2BC_ETA(H) \
	+ 1.2809002730136E-04 * SQ(REGION2_B2BC_ETA(H)) \
	)))

#endif


