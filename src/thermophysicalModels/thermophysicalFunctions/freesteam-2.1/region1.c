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
#include "region1.h"

static double gam(double pi, double tau);
static double gampi(double pi, double tau);
static double gampipi(double pi, double tau);
static double gamtau(double pi, double tau);
static double gamtautau(double pi, double tau);
static double gampitau(double pi, double tau);

#define REGION1_GPT_PSTAR 16.53e6 /* Pa */
#define REGION1_GPT_TSTAR 1386. /* K */

#define DEFINE_PITAU(P,T) \
	double pi = p / REGION1_GPT_PSTAR; \
	double tau = REGION1_GPT_TSTAR / T

#define R IAPWS97_R

#include <math.h>

double freesteam_region1_u_pT(double p, double T){
	DEFINE_PITAU(P,T);
	return (R * T) * (tau * gamtau(pi,tau) - pi * gampi(pi,tau));
}

double freesteam_region1_v_pT(double p, double T){
	DEFINE_PITAU(P,T);
	return (R * T / p) * pi * gampi(pi,tau);
}

double freesteam_region1_s_pT(double p, double T){
	DEFINE_PITAU(P,T);
	return R * (tau * gamtau(pi,tau) - gam(pi,tau));
}

double freesteam_region1_h_pT(double p, double T){
	DEFINE_PITAU(P,T);
	return R * T * (tau * gamtau(pi,tau));
}

double freesteam_region1_cp_pT(double p, double T){
	DEFINE_PITAU(P,T);
	return R * (-SQ(tau) * gamtautau(pi,tau));
}

double freesteam_region1_cv_pT(double p, double T){
	DEFINE_PITAU(P,T);
	return R * (-SQ(tau) * gamtautau(pi,tau) + SQ(gampi(pi,tau) -
		tau * gampitau(pi,tau)) / gampipi(pi,tau)
	);
}

double freesteam_region1_w_pT(double p, double T){
	DEFINE_PITAU(P,T);
	double gp = gampi(pi,tau);
	return sqrt(R * T * SQ(gp) / \
		(SQ(gp - tau*gampitau(pi,tau))/SQ(tau)/gamtautau(pi,tau) - gampipi(pi,tau))
	);
}

double freesteam_region1_g_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return R * T * gam(pi,tau);
}

double freesteam_region1_a_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return R * T * (gam(pi,tau) - gampi(pi,tau) * pi);
}


double freesteam_region1_alphav_pT(double p, double T){
	DEFINE_PITAU(P,T);
	return 1./T * (1. - tau*gampitau(pi,tau)/gampi(pi,tau));
}

double freesteam_region1_kappaT_pT(double p, double T){
	DEFINE_PITAU(P,T);
	return -1./p * pi*gampipi(pi,tau)/gampi(pi,tau);
}

//----------------------------------------------------------------
// REGION 1 G(p,T) EQUATIONS

typedef struct{
	int I, J;
	double n;
} IJNData;

const IJNData REGION1_GPT_DATA[] = {
	{0,	-2,	0.14632971213167E+00}
	,{0,	-1,	-0.84548187169114E+00}
	,{0,	0,	-0.37563603672040E+01}
	,{0,	1,	0.33855169168385E+01}
	,{0,	2,	-0.95791963387872E+00}
	,{0,	3,	0.15772038513228E+00}
	,{0,	4,	-0.16616417199501E-01}
	,{0,	5,	0.81214629983568E-03}
	,{1,	-9,	0.28319080123804E-03}
	,{1,	-7,	-0.60706301565874E-03}
	,{1,	-1,	-0.18990068218419E-01}
	,{1,	0,	-0.32529748770505E-01}
	,{1,	1,	-0.21841717175414E-01}
	,{1,	3,	-0.52838357969930E-04}
	,{2,	-3,	-0.47184321073267E-03}
	,{2,	0,	-0.30001780793026E-03}
	,{2,	1,	0.47661393906987E-04}
	,{2,	3,	-0.44141845330846E-05}
	,{2,	17,	-0.72694996297594E-15}
	,{3,	-4,	-0.31679644845054E-04}
	,{3,	0,	-0.28270797985312E-05}
	,{3,	6,	-0.85205128120103E-09}
	,{4,	-5,	-0.22425281908000E-05}
	,{4,	-2,	-0.65171222895601E-06}
	,{4,	10,	-0.14341729937924E-12}
	,{5,	-8,	-0.40516996860117E-06}
	,{8,	-11,	-0.12734301741641E-08}
	,{8,	-6,	-0.17424871230634E-09}
	,{21,	-29,	-0.68762131295531E-18}
	,{23,	-31,	0.14478307828521E-19}
	,{29,	-38,	0.26335781662795E-22}
	,{30,	-39,	-0.11947622640071E-22}
	,{31,	-40,	0.18228094581404E-23}
	,{32,	-41,	-0.93537087292458E-25}
};

const unsigned REGION1_GPT_MAX = sizeof(REGION1_GPT_DATA)/sizeof(IJNData);

#define REGION1_GPT_LOOP \
	double sum = 0; \
	const IJNData *d, *e = REGION1_GPT_DATA + REGION1_GPT_MAX; \
	for(d = REGION1_GPT_DATA; d < e; ++d)

double gam(double pi, double tau){
	REGION1_GPT_LOOP{
		sum += d->n * ipow(7.1 - pi,d->I) * ipow(tau - 1.222,d->J);
	}
	return sum;
}

double gampi(double pi, double tau){
	REGION1_GPT_LOOP{
		sum +=  -d->n * d->I * ipow(7.1 - pi,d->I -1) * ipow(tau - 1.222,d->J);
	}
	return sum;
}

double gampipi(double pi, double tau){
	REGION1_GPT_LOOP{
		sum += d->n * d->I * (d->I - 1) * ipow(7.1 - pi, d->I - 2) * ipow(tau - 1.222, d->J);
    }
	return sum;
}

double gamtau(double pi, double tau){
	REGION1_GPT_LOOP{
		sum += d->n * ipow(7.1 - pi, d->I) * d->J * ipow(tau - 1.222, d->J - 1);
	}
	return sum;
}

double gamtautau(double pi, double tau){
	REGION1_GPT_LOOP{
		sum += d->n * ipow(7.1 - pi, d->I) * d->J * (d->J - 1) * ipow(tau - 1.222, d->J - 2);
	}
	return sum;
}

double gampitau(double pi, double tau){
	REGION1_GPT_LOOP{
		sum += -d->n * d->I * ipow(7.1 - pi, d->I - 1) * d->J * ipow(tau - 1.222, d->J - 1);
	}
	return sum;
}

