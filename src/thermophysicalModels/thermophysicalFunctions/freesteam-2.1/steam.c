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
#include "steam.h"

#include <stdlib.h>
#include <stdio.h>

#include "region1.h"
#include "region2.h"
#include "region3.h"
#include "region4.h"
#include "b23.h"
#include "backwards.h"
#include "viscosity.h"
#include "thcond.h"

/* 'setter' functions for SteamState (forwards equations) */

SteamState freesteam_region1_set_pT(double p, double T){
	SteamState S;
	S.region = 1;
	S.R1.p = p;
	S.R1.T = T;
	/* FIXME add bounds check? */
	return S;
}

SteamState freesteam_region2_set_pT(double p, double T){
	SteamState S;
	S.region = 2;
	S.R2.p = p;
	S.R2.T = T;
	/* FIXME add bounds check? */
	return S;
}

SteamState freesteam_region3_set_rhoT(double rho, double T){
	SteamState S;
	S.region = 3;
	S.R3.rho = rho;
	S.R3.T = T;
	/* FIXME add bounds check? */
	return S;
}

SteamState freesteam_region4_set_Tx(double T, double x){
	SteamState S;
	S.region = 4;
	S.R4.T = T;
	S.R4.x = x;
	/* FIXME add bounds check? */
	return S;
}

int freesteam_fprint(FILE *f, SteamState S){
	int n = 0;
	n += fprintf(f, "region %d: ", S.region);
	switch(S.region){
		case 1:
			n += fprintf(f, "p = %f MPa, T = %f K\n", S.R1.p/1e6, S.R1.T);
			break;
		case 2:
			n += fprintf(f, "p = %f MPa, T = %f K\n", S.R2.p/1e6, S.R2.T);
			break;
		case 3:
			n += fprintf(f, "rho = %f kg/m³, T = %f K\n", S.R3.rho, S.R1.T);
			break;
		case 4:
			n += fprintf(f, "T = %f, x = %f\n", S.R4.T, S.R4.x);
			break;
	}
	return n;
}

/* 'getter' functions for SteamState */

int freesteam_region(SteamState S){
	return (int)S.region;
}

double freesteam_T(SteamState S){
	switch(S.region){
		case 1:
			return S.R1.T;
		case 2:
			return S.R2.T;
		case 3:
			return S.R3.T;
		case 4:
			return S.R4.T;
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_T\n");
			exit(1);
	}
}

double freesteam_p(SteamState S){
	switch(S.region){
		case 1:
			return S.R1.p;
		case 2:
			return S.R2.p;
		case 3:
			return freesteam_region3_p_rhoT(S.R3.rho, S.R3.T);
		case 4:
			return freesteam_region4_psat_T(S.R4.T);
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_p\n");
			exit(1);
	}
}


double freesteam_v(SteamState S){
	switch(S.region){
		case 1:
			return freesteam_region1_v_pT(S.R1.p,S.R1.T);
		case 2:
			return freesteam_region2_v_pT(S.R2.p,S.R2.T);
		case 3:
			return 1./S.R3.rho;
		case 4:
			return freesteam_region4_v_Tx(S.R4.T, S.R4.x);
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_v\n");
			exit(1);
	}
}

double freesteam_rho(SteamState S){
	switch(S.region){
		case 1:
			return 1./freesteam_region1_v_pT(S.R1.p,S.R1.T);
		case 2:
			return 1./freesteam_region2_v_pT(S.R2.p,S.R2.T);
		case 3:
			return S.R3.rho;
		case 4:
			return 1./freesteam_region4_v_Tx(S.R4.T, S.R4.x);
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_rho\n");
			exit(1);
	}
}


double freesteam_u(SteamState S){
	switch(S.region){
		case 1:
			return freesteam_region1_u_pT(S.R1.p, S.R1.T);
		case 2:
			return freesteam_region2_u_pT(S.R2.p, S.R2.T);
		case 3:
			return freesteam_region3_u_rhoT(S.R3.rho,S.R3.T);
		case 4:
			return freesteam_region4_u_Tx(S.R4.T, S.R4.x);
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_u\n");
			exit(1);
	}
}

double freesteam_h(SteamState S){
	switch(S.region){
		case 1:
			return freesteam_region1_h_pT(S.R1.p, S.R1.T);
		case 2:
			return freesteam_region2_h_pT(S.R2.p, S.R2.T);
		case 3:
			return freesteam_region3_h_rhoT(S.R3.rho,S.R3.T);
		case 4:
			return freesteam_region4_h_Tx(S.R4.T, S.R4.x);
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_h\n");
			exit(1);
	}
}


double freesteam_s(SteamState S){
	switch(S.region){
		case 1:
			return freesteam_region1_s_pT(S.R1.p, S.R1.T);
		case 2:
			return freesteam_region2_s_pT(S.R2.p, S.R2.T);
		case 3:
			return freesteam_region3_s_rhoT(S.R3.rho,S.R3.T);
		case 4:
			return freesteam_region4_s_Tx(S.R4.T, S.R4.x);
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_s\n");
			exit(1);
	}
}

double freesteam_cp(SteamState S){
	switch(S.region){
		case 1:
			return freesteam_region1_cp_pT(S.R1.p, S.R1.T);
		case 2:
			return freesteam_region2_cp_pT(S.R2.p, S.R2.T);
		case 3:
			return freesteam_region3_cp_rhoT(S.R3.rho,S.R3.T);
		case 4:
			return freesteam_region4_cp_Tx(S.R4.T, S.R4.x);
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_cp\n");
			exit(1);
	}
}

double freesteam_cv(SteamState S){
	switch(S.region){
		case 1:
			return freesteam_region1_cv_pT(S.R1.p, S.R1.T);
		case 2:
			return freesteam_region2_cv_pT(S.R2.p, S.R2.T);
		case 3:
			return freesteam_region3_cv_rhoT(S.R3.rho,S.R3.T);
		case 4:
			return freesteam_region4_cv_Tx(S.R4.T, S.R4.x);
		default:
			fprintf(stderr,"ERROR: invalid region in freesteam_cv\n");
			exit(1);
	}
}

double freesteam_w(SteamState S){
	switch(S.region){
		case 1:
			return freesteam_region1_w_pT(S.R1.p, S.R1.T);
		case 2:
			return freesteam_region2_w_pT(S.R2.p, S.R2.T);
		case 3:
			return freesteam_region3_w_rhoT(S.R3.rho,S.R3.T);
#if 0
		case 4:
			return freesteam_region4_w_Tx(S.R4.T, S.R4.x);
#endif
		default:
			fprintf(stderr,"ERROR: invalid region '%d' in freesteam_w\n",S.region);
			exit(1);
	}
}

double freesteam_x(SteamState S){
	switch(S.region){
		case 1:
			return 0.;
		case 2:
			return 1.;
		case 3:
			if(S.R3.rho > IAPWS97_RHOCRIT)return 0.;
			return 1.;
		case 4:
			return S.R4.x;
		default:
			fprintf(stderr,"ERROR: invalid region '%d' in freesteam_x\n",S.region);
			exit(1);
	}
}

double freesteam_mu(SteamState S){
	static char warned = 0;
	switch(S.region){
		case 1:
            return freesteam_mu_rhoT(1./freesteam_region1_v_pT(S.R1.p,S.R1.T), S.R1.T);
		case 2:
            return freesteam_mu_rhoT(1./freesteam_region2_v_pT(S.R2.p,S.R2.T), S.R2.T);
		case 3:
            return freesteam_mu_rhoT(S.R3.rho, S.R3.T);
		case 4:
			if(!warned){
				fprintf(stderr,"WARNING: viscosity evaluation in region 4 is poorly defined! (this warning is only shown once)\n");
				warned = 1;
			}
            return freesteam_mu_rhoT(1./freesteam_region4_v_Tx(S.R4.T, S.R4.x), S.R4.T);
		default:
			fprintf(stderr,"ERROR: invalid region '%d' in freesteam_mu\n",S.region);
			exit(1);
	}
}

double freesteam_k(SteamState S){
	static char warned = 0;
	switch(S.region){
		case 1:
            return freesteam_k_rhoT(1./freesteam_region1_v_pT(S.R1.p,S.R1.T), S.R1.T);
		case 2:
            return freesteam_k_rhoT(1./freesteam_region2_v_pT(S.R2.p,S.R2.T), S.R2.T);
		case 3:
            return freesteam_k_rhoT(S.R3.rho, S.R3.T);
		case 4:
			if(!warned){
				fprintf(stderr,"WARNING: thermal conductivity evaluation in region 4 is poorly defined! (this warning is only shown once)\n");
				warned = 1;
			}
            return freesteam_k_rhoT(1./freesteam_region4_v_Tx(S.R4.T, S.R4.x), S.R4.T);
		default:
			fprintf(stderr,"ERROR: invalid region '%d' in freesteam_k\n",S.region);
			exit(1);
	}
}
