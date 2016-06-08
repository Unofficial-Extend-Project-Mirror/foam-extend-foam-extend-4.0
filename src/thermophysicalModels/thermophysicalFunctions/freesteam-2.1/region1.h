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

#ifndef FREESTEAM_REGION1_H
#define FREESTEAM_REGION1_H

#include "common.h"

FREESTEAM_DLL double freesteam_region1_u_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_v_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_s_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_h_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_cp_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_cv_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_w_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_a_pT(double p, double T);
FREESTEAM_DLL double freesteam_region1_g_pT(double p, double T);

/* used in calculations of derivatives, see derivs.c */
double freesteam_region1_alphav_pT(double p, double T);
double freesteam_region1_kappaT_pT(double p, double T);

#define REGION1_TMAX 623.15 /* K */

#endif

