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

/* Appendix B: Recommended Interpolating equation for Industrial Use */
/* see http://www.iapws.org/relguide/thcond.pdf */

#define FREESTEAM_BUILDING_LIB
#include "thcond.h"

#include <math.h>

#define THCOND_TSTAR 647.26
#define THCOND_RHOSTAR 317.7
#define THCOND_KSTAR 1.0

#define THCOND_b0 -0.397070
#define THCOND_b1 0.400302
#define THCOND_b2 1.060000
#define THCOND_B1 -0.171587
#define THCOND_B2 2.392190

#define THCOND_C1 0.642857
#define THCOND_C2 -4.11717
#define THCOND_C3 -6.17937
#define THCOND_C4 0.00308976
#define THCOND_C5 0.0822994
#define THCOND_C6 10.0932

#define THCOND_d1 0.0701309
#define THCOND_d2 0.0118520
#define THCOND_d3 0.00169937
#define THCOND_d4 -1.0200

/* freesteam code */
double freesteam_k_rhoT(double rho, double T){

#define THCOND_a_COUNT 4
	const double THCOND_a[THCOND_a_COUNT] = {
		0.0102811
		,0.0299621
		,0.0156146
		,-0.00422464
	};

	double Tbar = T / THCOND_TSTAR;
	double rhobar = rho / THCOND_RHOSTAR;

	/* fast implementation... minimised calls to 'pow' routine... */

	double Troot = sqrt(Tbar);
	double Tpow = Troot;
	double lam = 0;

	int k;
	for(k = 0; k < THCOND_a_COUNT; ++k) {
		lam += THCOND_a[k] * Tpow;
		Tpow *= Tbar;
	}

	lam += THCOND_b0 + THCOND_b1 * rhobar + THCOND_b2 * exp(THCOND_B1 * SQ(rhobar + THCOND_B2));
	
	double DTbar = fabs(Tbar - 1) + THCOND_C4;
	double DTbarpow = pow(DTbar, 3./5);
	double Q = 2. + THCOND_C5 / DTbarpow;

	double S;
	if(Tbar >= 1){
		S = 1. / DTbar;
	}else{
		S = THCOND_C6 / DTbarpow;
	}

	double rhobar18 = pow(rhobar, 1.8);
	double rhobarQ = pow(rhobar, Q);

	lam += 
		(THCOND_d1 / ipow(Tbar,10) + THCOND_d2) * rhobar18 * 
			exp(THCOND_C1 * (1 - rhobar * rhobar18))
		+ THCOND_d3 * S * rhobarQ *
			exp((Q/(1+Q))*(1 - rhobar*rhobarQ))
		+ THCOND_d4 *
			exp(THCOND_C2 * ipow(Troot,3) + THCOND_C3 / ipow(rhobar,5));

	return THCOND_KSTAR * lam;
}

