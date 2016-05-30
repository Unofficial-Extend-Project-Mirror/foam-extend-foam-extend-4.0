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
#include "region2.h"

#define GAM0(PI,TAU) gam0(PI,TAU)
#define GAM0PI(PI,TAU) (1./PI)
#define GAM0PIPI(PI,TAU) (-1./SQ(PI))
#define GAM0TAU(PI,TAU) gam0tau(TAU)
#define GAM0TAUTAU(PI,TAU) gam0tautau(TAU)
#define GAM0PITAU(PI,TAU) (0)

#define gam(PI,TAU) (GAM0(PI,TAU) + gamr(PI,TAU))
#define gampi(PI,TAU) (GAM0PI(PI,TAU) + gamrpi(PI,TAU))
#define gampipi(PI,TAU) (GAM0PIPI(PI,TAU) + gamrpipi(PI,TAU))
#define gamtau(PI,TAU) (GAM0TAU(PI,TAU) + gamrtau(PI,TAU))
#define gamtautau(PI,TAU) (GAM0TAUTAU(PI,TAU) + gamrtautau(PI,TAU))
#define gampitau(PI,TAU) (GAM0PITAU(PI,TAU) + gamrpitau(PI,TAU))

static double gamr(double pi, double tau);
static double gamrpi(double pi, double tau);
static double gamrpipi(double pi, double tau);
static double gamrtau(double pi, double tau);
static double gamrtautau(double pi, double tau);
static double gamrpitau(double pi, double tau);

static double gam0(double pi, double tau);
static double gam0tau(double tau);
static double gam0tautau(double tau);

#include <math.h>
#include "common.h"

#define R IAPWS97_R

#define REGION2_GPT_PSTAR 1.e6 /* Pa */
#define REGION2_GPT_TSTAR 540. /* K */

#define DEFINE_PITAU(P,T) \
	double pi = p / REGION2_GPT_PSTAR; \
	double tau = REGION2_GPT_TSTAR / T


double freesteam_region2_v_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return (R * T / p) * pi * gampi(pi,tau);
}

double freesteam_region2_u_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return (R * T) * (tau * gamtau(pi,tau) - pi * gampi(pi,tau));
}

double freesteam_region2_s_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return R * (tau * gamtau(pi,tau) - gam(pi,tau));
}

double freesteam_region2_h_pT(double p, double T){
	//fprintf(stderr,"%s: p = %f, T = %f\n",__func__,p,T);
	DEFINE_PITAU(p,T);
	return R * T * (tau * gamtau(pi,tau));
}

double freesteam_region2_cp_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return R * (-SQ(tau) * gamtautau(pi,tau));
}

double freesteam_region2_cv_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return R * (-SQ(tau) * gamtautau(pi,tau) + SQ(gampi(pi,tau) -
		tau * gampitau(pi,tau)) / gampipi(pi,tau)
	);
}

double freesteam_region2_w_pT(double p, double T){
	DEFINE_PITAU(p,T);
	double gp = gamrpi(pi,tau);
	return sqrt(R * T * (1. + 2.*pi*gp+SQ(pi*gp))/
		((1. - SQ(pi)*gamrpipi(pi,tau))	+ SQ(1. + pi*gp - tau*pi*gamrpitau(pi,tau))/SQ(tau)/gamtautau(pi,tau))
	);
}

double freesteam_region2_g_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return R * T * gam(pi,tau);
}	

double freesteam_region2_a_pT(double p, double T){
	DEFINE_PITAU(p,T);
	return R * T * (gam(pi,tau) - gampi(pi,tau) * pi);
}

double freesteam_region2_alphav_pT(double p, double T){
	DEFINE_PITAU(p,T);
	double pigamrpi = pi*gamrpi(pi,tau);
	double alphav = 1./T * (1. + pigamrpi - tau*pi*gamrpitau(pi,tau))/(1. + pigamrpi);
	//fprintf(stderr,"α_v = %g\n",alphav);
	return alphav;
}

double freesteam_region2_kappaT_pT(double p, double T){
	DEFINE_PITAU(p,T);
	double kappaT = 1./p * (1.-SQ(pi)*gamrpipi(pi,tau)) / (1.+pi*gamrpi(pi,tau));
	//fprintf(stderr,"κ_T = %g\n",kappaT);
	return kappaT;
}

/*------------------------------------------------------------------------------
  REGION 2 IDEAL PART - GAM0(PI,TAU)
*/

typedef struct{
	int J;
	double n;
} JNData;

const JNData REGION2_GPT_IDEAL_DATA[] = {
	{0,	-0.96927686500217E+01}
	,{1,	0.10086655968018E+02}
	,{-5,	-0.56087911283020E-02}
	,{-4,	0.71452738081455E-01}
	,{-3,	-0.40710498223928E+00}
	,{-2,	0.14240819171444E+01}
	,{-1,	-0.43839511319450E+01}
	,{2,	-0.28408632460772E+00}
	,{3,	0.21268463753307E-01}
};

const unsigned REGION2_GPT_IDEAL_MAX = sizeof(REGION2_GPT_IDEAL_DATA)/sizeof(JNData);

#define REGION2_GPT_IDEAL_LOOP \
	double sum = 0; \
	const JNData *d, *e = REGION2_GPT_IDEAL_DATA + REGION2_GPT_IDEAL_MAX; \
	for(d = REGION2_GPT_IDEAL_DATA; d < e; ++d)

double gam0(double pi, double tau){
	REGION2_GPT_IDEAL_LOOP{
		sum += d->n * ipow(tau, d->J);
	}
	return log(pi) + sum;
}

double gam0tau(double tau){
	REGION2_GPT_IDEAL_LOOP{
		sum += d->n * d->J * ipow(tau, d->J - 1);
	}
	return sum;
}

double gam0tautau(double tau){
	REGION2_GPT_IDEAL_LOOP{
		sum += d->n * d->J * (d->J - 1) * ipow(tau, d->J - 2);
	}
	return sum;
}

/*------------------------------------------------------------------------------
  REGION 2 RESIDUAL PART - GAMR(PI,TAU)
*/

typedef struct{
	int I, J;
	double n;
} IJNData;

const IJNData REGION2_GPT_RESID_DATA[] = {
	{1,	0,	-0.17731742473213E-02}
	,{1,	1,	-0.17834862292358E-01}
	,{1,	2,	-0.45996013696365E-01}
	,{1,	3,	-0.57581259083432E-01}
	,{1,	6,	-0.50325278727930E-01}
	,{2,	1,	-0.33032641670203E-04}
	,{2,	2,	-0.18948987516315E-03}
	,{2,	4,	-0.39392777243355E-02}
	,{2,	7,	-0.43797295650573E-01}
	,{2,	36,	-0.26674547914087E-04}
	,{3,	0,	0.20481737692309E-07}
	,{3,	1,	0.43870667284435E-06}
	,{3,	3,	-0.32277677238570E-04}
	,{3,	6,	-0.15033924542148E-02}
	,{3,	35,	-0.40668253562649E-01}
	,{4,	1,	-0.78847309559367E-09}
	,{4,	2,	0.12790717852285E-07}
	,{4,	3,	0.48225372718507E-06}
	,{5,	7,	0.22922076337661E-05}
	,{6,	3,	-0.16714766451061E-10}
	,{6,	16,	-0.21171472321355E-02}
	,{6,	35,	-0.23895741934104E+02}
	,{7,	0,	-0.59059564324270E-17}
	,{7,	11,	-0.12621808899101E-05}
	,{7,	25,	-0.38946842435739E-01}
	,{8,	8,	0.11256211360459E-10}
	,{8,	36,	-0.82311340897998E+01}
	,{9,	13,	0.19809712802088E-07}
	,{10,	4,	0.10406965210174E-18}
	,{10,	10,	-0.10234747095929E-12}
	,{10,	14,	-0.10018179379511E-08}
	,{16,	29,	-0.80882908646985E-10}
	,{16,	50,	0.10693031879409E+00}
	,{18,	57,	-0.33662250574171E+00}
	,{20,	20,	0.89185845355421E-24}
	,{20,	35,	0.30629316876232E-12}
	,{20,	48,	-0.42002467698208E-05}
	,{21,	21,	-0.59056029685639E-25}
	,{22,	53,	0.37826947613457E-05}
	,{23,	39,	-0.12768608934681E-14}
	,{24,	26,	0.73087610595061E-28}
	,{24,	40,	0.55414715350778E-16}
	,{24,	58,	-0.94369707241210E-06}
};

const unsigned REGION2_GPT_RESID_MAX = sizeof(REGION2_GPT_RESID_DATA)/sizeof(IJNData);

#define REGION2_GPT_RESID_LOOP \
	double sum = 0; \
	const IJNData *d, *e = REGION2_GPT_RESID_DATA + REGION2_GPT_RESID_MAX; \
	for(d = REGION2_GPT_RESID_DATA; d < e; ++d)

double gamr(double pi, double tau){
	REGION2_GPT_RESID_LOOP{
		sum += d->n * ipow(pi,d->I) * ipow(tau - 0.5,d->J);
	}
	return sum;
}

double gamrpi(double pi, double tau){
	REGION2_GPT_RESID_LOOP{
		sum +=  d->n * d->I * ipow(pi,d->I -1) * ipow(tau - 0.5,d->J);
	}
	return sum;
}

double gamrpipi(double pi, double tau){
	REGION2_GPT_RESID_LOOP{
		sum += d->n * d->I * (d->I - 1) * ipow(pi, d->I - 2) * ipow(tau - 0.5, d->J);
    }
	return sum;
}

double gamrtau(double pi, double tau){
	REGION2_GPT_RESID_LOOP{
		sum += d->n * ipow(pi, d->I) * d->J * ipow(tau - 0.5, d->J - 1);
	}
	return sum;
}

double gamrtautau(double pi, double tau){
	REGION2_GPT_RESID_LOOP{
		sum += d->n * ipow(pi, d->I) * d->J * (d->J - 1) * ipow(tau - 0.5, d->J - 2);
	}
	return sum;
}

double gamrpitau(double pi, double tau){
	REGION2_GPT_RESID_LOOP{
		sum += d->n * d->I * ipow(pi, d->I - 1) * d->J * ipow(tau - 0.5, d->J - 1);
	}
	return sum;
}


