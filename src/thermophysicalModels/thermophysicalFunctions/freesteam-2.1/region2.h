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

#ifndef FREESTEAM_REGION2_H
#define FREESTEAM_REGION2_H

#include "common.h"

#define REGION2_TMAX 1073.15

FREESTEAM_DLL double freesteam_region2_v_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_u_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_s_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_h_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_cp_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_cv_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_w_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_a_pT(double p, double T);
FREESTEAM_DLL double freesteam_region2_g_pT(double p, double T);

/* used in calculations of derivatives, see derivs.c */
double freesteam_region2_alphav_pT(double p, double T);
double freesteam_region2_kappaT_pT(double p, double T);

#endif

