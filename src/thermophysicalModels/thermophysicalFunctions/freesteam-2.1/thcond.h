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

#ifndef FREESTEAM_THCOND_H
#define FREESTEAM_THCOND_H

#include "common.h"

/// Conductivity [W/m.K]
/**
	Returns the thermal conductivity of water/steam.
	@see http://www.iapws.org/relguide/thcond.pdf

	Range of validity is entire regions 1,2,3. The correlation is not
	really applicable in region 4, but will give 'sane' results there.

	@return Thermal conductivity [W/m.K]
*/
FREESTEAM_DLL double freesteam_k_rhoT(double rho, double T);

#endif

