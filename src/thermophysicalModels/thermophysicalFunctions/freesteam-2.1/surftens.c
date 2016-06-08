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
*/
#define FREESTEAM_BUILDING_LIB
#include "surftens.h"
#include <math.h>

/**
	Calculate surface tension between water phase and vapour phase

	@param T temperature (Kelvin)
	@return surface tension (N/m)

	The correlation is the IAPWS Release on Surface Tension of Ordinary Water
	Substance, September 1994.
	@since 0.7
*/
double freesteam_surftens_T(double T){
	double tau = 1. - T / IAPWS97_TCRIT;
	const double B = 235.8e-3; /* converted to N/m */
	const double b = -0.625;
	const double mu = 1.256;

	return B * pow(tau,mu) * (1 + b * tau);
}

