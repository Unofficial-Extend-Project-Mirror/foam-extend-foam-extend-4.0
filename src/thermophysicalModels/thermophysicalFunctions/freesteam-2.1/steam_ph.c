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
#include "steam_ph.h"

#include "region1.h"
#include "region2.h"
#include "region3.h"
#include "region4.h"
#include "b23.h"
#include "backwards.h"
#include "solver2.h"
#include "zeroin.h"

#include "common.h"

#include <stdio.h>
#include <stdlib.h>

int freesteam_bounds_ph(double p, double h, int verbose){

#define BOUND_WARN(MSG) \
	if(verbose){\
		fprintf(stderr,"%s (%s:%d): WARNING " MSG " (p = %g MPa, h = %g kJ/kg)\n"\
		,__func__,__FILE__,__LINE__,p/1e6,h/1e3);\
	}

	if(p <= 0){
		BOUND_WARN("p <= 0");
		return 1;
	}
	if(p > IAPWS97_PMAX){
		BOUND_WARN("p > PMAX");
		return 2;
	}

	double hmax = freesteam_region2_h_pT(p,REGION2_TMAX);
	if(h>hmax){
		BOUND_WARN("h > hmax");
		return 3;
	}
	double hmin = freesteam_region1_h_pT(p,IAPWS97_TMIN);
	if(h < hmin){
		BOUND_WARN("h < hmin");
		return 4;
	}
	return 0;
#undef BOUND_WARN
}

int freesteam_region_ph(double p, double h){
	//fprintf(stderr,"freesteam_set_ph(p = %f, h = %f)\n",p,h);

	double p13 = 0;
	p13 = freesteam_region4_psat_T(REGION1_TMAX);

	//fprintf(stderr,"p13 = %lf MPa\n",p13/1.e6);
	//fprintf(stderr,"check: %f\n",freesteam_region4_psat_T(REGION1_TMAX));

	if(p <= p13){
		double Tsat = freesteam_region4_Tsat_p(p);
		double hf = freesteam_region1_h_pT(p,Tsat);
		if(h<hf){
			return 1;
		}
		double hg = freesteam_region2_h_pT(p,Tsat);
		if(h>hg){
			return 2;
		}
		/* this is the low-pressure portion of region 4 */
		return 4;
	}

	double h13 = freesteam_region1_h_pT(p,REGION1_TMAX);
	if(h <= h13){
		return 1;
	}

	double T23 = freesteam_b23_T_p(p);
	//fprintf(stderr,"p = %f MPa --> T23(p) = %f K (%f °C)\n",p/1e6,T23,T23-273.15);
	double h23 = freesteam_region2_h_pT(p,T23);
	if(h >= h23){
		return 2;
	}

	double psat = freesteam_region3_psat_h(h);
	if(p > psat){
		return 3;
	}

	/* high-pressure portion of region 4 */
	return 4;
}


typedef struct SolvePHData_struct{
	double p, h;
} SolvePHData;

#define D ((SolvePHData *)user_data)
static ZeroInSubjectFunction ph_region2_fn;
double ph_region2_fn(double T, void *user_data){
	return D->h - freesteam_region2_h_pT(D->p, T);
}
#undef D


SteamState freesteam_set_ph(double p, double h){
	SteamState S;
	S.region = (char)freesteam_region_ph(p,h);
	//int status;
	switch(S.region){
		case 1:
			S.R1.p = p;
			S.R1.T = freesteam_region1_T_ph(p, h);
#if 0
			S = freesteam_solver2_region1('p','h', p, h, S, &status);
			if(status){
				fprintf(stderr,"%s: WARNING: Failed to converge in region 1\n",__func__);
			}
#endif
			return S;
		case 2:
			S.R2.p = p;
			S.R2.T = freesteam_region2_T_ph(p, h);
			{
				double lb = S.R2.T * 0.999;
				double ub = S.R2.T * 1.001;
				double tol = 1e-9; /* ??? */
				double sol, err;
				SolvePHData D = {p, h};
				zeroin_solve(&ph_region2_fn, &D, lb, ub, tol, &sol, &err);
				S.R2.T = sol;
			}
#if 0
			/* solver2 is not working in this region, for some reason. */
			S = freesteam_solver2_region2('p','h', p, h, S, &status);
			if(status){
				fprintf(stderr,"%s: WARNING: Failed to converge in region 2\n");
			}
#endif
			return S;
		case 3:
			S.R3.rho = 1./freesteam_region3_v_ph(p, h);
			S.R3.T = freesteam_region3_T_ph(p, h);
#if 0
			/* FIXME: this code doesn't work, see pTdiagram.h for example */
			/* by not iterating here, the relative error is increased from
			about 1e-12 to 1e-4 */
			S = freesteam_solver2_region3('p','h', p, h, S, &status);
			if(status){
				fprintf(stderr,"%s: WARNING: Failed to converge in region 3\n",__func__);
			}
#endif
			return S;
		case 4:
			S.R4.T = freesteam_region4_Tsat_p(p);
			//fprintf(stderr,"%s: region 4, Tsat = %g\n",__func__,S.R4.T);
			double hf, hg;
			if(S.R4.T <= REGION1_TMAX){
				hf = freesteam_region1_h_pT(p,S.R4.T);
				hg = freesteam_region2_h_pT(p,S.R4.T);
			}else{
				/* TODO iteratively improve estimate of T */
				double rhof = freesteam_region4_rhof_T(S.R4.T);
				double rhog = freesteam_region4_rhog_T(S.R4.T);
				/* FIXME iteratively improve these estimates of rhof, rhog */
				hf = freesteam_region3_h_rhoT(rhof,S.R4.T);
				hg = freesteam_region3_h_rhoT(rhog,S.R4.T);
			}
			S.R4.x = (h - hf)/(hg - hf);
			return S;
		default:
			fprintf(stderr,"%s: invalid region %d\n",__func__,S.region);
			return S;
	}
}




