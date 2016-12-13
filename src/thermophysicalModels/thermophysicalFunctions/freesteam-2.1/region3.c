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
#include "region3.h"

const double REGION3_ARHOT_TSTAR = 647.096 /* K */;
const double REGION3_ARHOT_RHOSTAR = 322. /* K */;

#define DEFINE_DELTAU(RHO,T) \
	double del = rho / REGION3_ARHOT_RHOSTAR; \
	double tau = REGION3_ARHOT_TSTAR / T

#define R 461.526

static double phi(double del, double tau);
static double phidel(double del, double tau);
static double phideldel(double del, double tau);
static double phitau(double del, double tau);
static double phitautau(double del, double tau);
static double phideltau(double del, double tau);

#include <math.h>

double freesteam_region3_p_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return rho * R * T * del * phidel(del,tau);
}

double freesteam_region3_u_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return R * T * tau * phitau(del,tau);
}

double freesteam_region3_s_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return R * (tau * phitau(del,tau) - phi(del,tau));
}

double freesteam_region3_h_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return R * T * (tau * phitau(del,tau) + del * phidel(del,tau));
}

double freesteam_region3_cp_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return R * (
		-SQ(tau) * phitautau(del,tau)
		+ (
			ipow (del * phidel(del,tau) - del * tau * phideltau(del,tau), 2)
			/ (2 * del * phidel(del,tau) + SQ(del) * phideldel(del,tau))
		)
	);
}

double freesteam_region3_cv_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return R * (-SQ(tau) * phitautau(del,tau));
}

double freesteam_region3_w_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return sqrt(R * T * (
		2 * del * phidel(del,tau) + SQ(del) * phideldel(del,tau)
		- (
			ipow (del * phidel(del,tau) - del * tau * phideltau(del,tau), 2)
			/ (SQ(tau) * phitautau(del,tau))
		)
	));
}

double freesteam_region3_alphap_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return 1./T * (1. - tau*phideltau(del,tau)/phidel(del,tau));
}

double freesteam_region3_betap_rhoT(double rho, double T){
	DEFINE_DELTAU(rho,T);
	return rho*(2. + del * phideldel(del,tau)/phidel(del,tau));
}

/*----------------------------------------------------------------------------*/

typedef struct{
	int I, J;
	double n;
} IJNData;

const double REGION3_N1 = 0.10658070028513E+01;

const IJNData REGION3_ARHOT_DATA[] = {
	{0,	0,	-0.15732845290239E+02}
	,{0,	1,	0.20944396974307E+02}
	,{0,	2,	-0.76867707878716E+01}
	,{0,	7,	0.26185947787954E+01}
	,{0,	10,	-0.28080781148620E+01}
	,{0,	12,	0.12053369696517E+01}
	,{0,	23,	-0.84566812812502E-02}
	,{1,	2,	-0.12654315477714E+01}
	,{1,	6,	-0.11524407806681E+01}
	,{1,	15,	0.88521043984318E+00}
	,{1,	17,	-0.64207765181607E+00}
	,{2,	0,	0.38493460186671E+00}
	,{2,	2,	-0.85214708824206E+00}
	,{2,	6,	0.48972281541877E+01}
	,{2,	7,	-0.30502617256965E+01}
	,{2,	22,	0.39420536879154E-01}
	,{2,	26,	0.12558408424308E+00}
	,{3,	0,	-0.27999329698710E+00}
	,{3,	2,	0.13899799569460E+01}
	,{3,	4,	-0.20189915023570E+01}
	,{3,	16,	-0.82147637173963E-02}
	,{3,	26,	-0.47596035734923E+00}
	,{4,	0,	0.43984074473500E-01}
	,{4,	2,	-0.44476435428739E+00}
	,{4,	4,	0.90572070719733E+00}
	,{4,	26,	0.70522450087967E+00}
	,{5,	1,	0.10770512626332E+00}
	,{5,	3,	-0.32913623258954E+00}
	,{5,	26,	-0.50871062041158E+00}
	,{6,	0,	-0.22175400873096E-01}
	,{6,	2,	0.94260751665092E-01}
	,{6,	26,	0.16436278447961E+00}
	,{7,	2,	-0.13503372241348E-01}
	,{8,	26,	-0.14834345352472E-01}
	,{9,	2,	0.57922953628084E-03}
	,{9,	26,	0.32308904703711E-02}
	,{10,	0,	0.80964802996215E-04}
	,{10,	1,	-0.16557679795037E-03}
	,{11,	26,	-0.44923899061815E-04}
};

const unsigned REGION3_ARHOT_MAX = sizeof(REGION3_ARHOT_DATA)/sizeof(IJNData);

#define REGION3_ARHOT_LOOP \
	double sum = 0; \
	const IJNData *d, *e = REGION3_ARHOT_DATA + REGION3_ARHOT_MAX; \
	for(d = REGION3_ARHOT_DATA; d < e; ++d)

double phi(double del, double tau){
	REGION3_ARHOT_LOOP{
		sum += d->n * ipow(del, d->I) * ipow(tau, d->J);
	}
	return sum + REGION3_N1 * log(del);
}

double phidel(double del, double tau){
	REGION3_ARHOT_LOOP{
		sum += +d->n * d->I * ipow(del, d->I - 1) * ipow(tau, d->J);
	}
	return sum + REGION3_N1 / del;
}

double phideldel(double del, double tau){
	REGION3_ARHOT_LOOP{
		sum += d->n * d->I * (d->I - 1) * ipow(del, d->I - 2) * ipow(tau, d->J);
	}
	return sum - REGION3_N1 / SQ(del) ;
}

double phitau(double del, double tau){
	REGION3_ARHOT_LOOP{
		sum += d->n * ipow(del, d->I) * d->J * ipow(tau, d->J - 1);
	}
	return sum;
}

double phitautau(double del, double tau){
	REGION3_ARHOT_LOOP{
		sum += d->n * ipow(del, d->I) * d->J * (d->J - 1) * ipow(tau, d->J - 2);
	}
	return sum;
}

double phideltau(double del, double tau){
	REGION3_ARHOT_LOOP{
		sum += d->n * d->I * ipow(del, d->I - 1) * d->J * ipow(tau, d->J - 1);
	}
	return sum;
}

