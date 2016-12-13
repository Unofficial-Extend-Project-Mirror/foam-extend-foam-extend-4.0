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
#ifndef FREESTEAM_BACKWARDS_H
#define FREESTEAM_BACKWARDS_H

#include "common.h"

FREESTEAM_DLL double freesteam_region1_T_ph(double p, double h);
FREESTEAM_DLL double freesteam_region2_T_ph(double p, double h);
FREESTEAM_DLL double freesteam_region3_T_ph(double p, double h);
FREESTEAM_DLL double freesteam_region3_v_ph(double p, double h);
FREESTEAM_DLL double freesteam_region3_psat_h(double h);
FREESTEAM_DLL double freesteam_region3_psat_s(double s);

FREESTEAM_DLL double freesteam_region3_T_ps(double p, double h);
FREESTEAM_DLL double freesteam_region3_v_ps(double p, double h);

#endif

