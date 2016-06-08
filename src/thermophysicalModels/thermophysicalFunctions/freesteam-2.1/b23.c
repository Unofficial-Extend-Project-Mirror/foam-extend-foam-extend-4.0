/*
freesteam - IAPWS-IF97 steam tables library
Copyright (C) 2004-2005  John Pye

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
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.*/

#define FREESTEAM_BUILDING_LIB
#include "b23.h"

#include <math.h>

const double B23_N[6] = { 0, 0.34805185628969E+03, -0.11671859879975E+01
	, 0.10192970039326E-02, 0.57254459862746E+03, 0.13918839778870E+02
};

const double B23_PSTAR = 1e6;

double freesteam_b23_p_T(double T){
#define theta T
	return (B23_N[1] + B23_N[2] * theta + B23_N[3] * SQ(theta)) * B23_PSTAR;
#undef theta
}

double freesteam_b23_T_p(double p){
	double pi = p / B23_PSTAR;
	return B23_N[4] + sqrt((pi - B23_N[5])/B23_N[3]) /* * 1{K} */;
}

