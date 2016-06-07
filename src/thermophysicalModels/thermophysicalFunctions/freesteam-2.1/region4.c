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
#include "region4.h"

#include "region1.h"
#include "region3.h"
#include "region2.h"

#include <math.h>
#include <stdlib.h>

const double REGION4_N[11] = { 0, 0.11670521452767E+04, -0.72421316703206E+06
	, -0.17073846940092E+02, 0.12020824702470E+05, -0.32325550322333E+07
	, 0.14915108613530E+02, -0.48232657361591E+04, 0.40511340542057E+06
	, -0.23855557567849E+00, 0.65017534844798E+03
};

#define REGION4_PSTAR 1e6 /* Pa */
#define REGION4_TSTAR 1 /* K */

/*------------------------------------------------------------------------------
  REGION 4 SATURATION CURVE psat(T)
*/

double freesteam_region4_psat_T(double T){

	//fprintf(stderr,"freesteam_region4_psat_T(T = %f)\n", T );

	double ups = T/REGION4_TSTAR + REGION4_N[9] / (T/REGION4_TSTAR - REGION4_N[10]);
	double A = SQ(ups) + REGION4_N[1] * ups + REGION4_N[2];
	double B = REGION4_N[3] * SQ(ups) + REGION4_N[4] * ups + REGION4_N[5];
	double C = REGION4_N[6] * SQ(ups) + REGION4_N[7] * ups + REGION4_N[8];

	double expr = 2. * C / (- B + sqrt(SQ(B) - 4. * A * C));
	double psat = SQ(SQ(expr)) * REGION4_PSTAR;

	/* fprintf(stderr,"freesteam_region4_psat_T = %f MPa\n", psat/1e6);*/
	return psat;
}

/*------------------------------------------------------------------------------
  REGION 4 SATURATION CURVE Tsat(p) (BACKWARDS EQUATION)
*/

double freesteam_region4_Tsat_p(double p){

	IAPWS97_APPROXIMATE;

	double beta = pow(p/REGION4_PSTAR, 0.25);
	double E = SQ(beta) + REGION4_N[3] * beta + REGION4_N[6];
	double F = REGION4_N[1] * SQ(beta) + REGION4_N[4] * beta + REGION4_N[7];
	double G = REGION4_N[2] * SQ(beta) + REGION4_N[5] * beta + REGION4_N[8];
	double D = 2. * G / (-F - sqrt(SQ(F) - 4. * E * G));

	double theta = 0.5 * (REGION4_N[10] + D - sqrt(SQ(REGION4_N[10] + D) - 4.0 * (REGION4_N[9] + REGION4_N[10] * D)));

	/* FIXME iterative improve this estimate? is it necessary? */

	return theta /* * REGION4_TSTAR = 1 {K} */;
}

/*------------------------------------------------------------------------------
  REGION 4 DENSITIES rhof(T), rhog(T) (SUPPLEMENTARY EQUATIONS)
*/



/**
	Coefficients for getSatDensWater_T
*/
const double REGION4_B[7]
	= { 0, 1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450E+05 };

/**
	Coefficients for getSatDensSteam_T
*/
const double REGION4_C[7]
	= { 0, -2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063 };


double freesteam_region4_rhof_T(double T){

	IAPWS97_APPROXIMATE;

	double tau = 1 - T / IAPWS97_TCRIT;

	double tau_1_3 = pow(tau,1./3);

	double tau_2_3 = SQ(tau_1_3);
	double tau_5_3 = tau * tau_2_3;
	double tau_16_3 = SQ(tau_5_3) * tau_5_3 * tau_1_3;
	double tau_43_3 = SQ(tau_16_3) * SQ(tau_5_3) * tau_1_3;
	double tau_110_3 = SQ(tau_43_3) * tau_16_3 * tau_5_3 * tau;

	double delta = 1
		+ REGION4_B[1]*tau_1_3
		+ REGION4_B[2]*tau_2_3
		+ REGION4_B[3]*tau_5_3
		+ REGION4_B[4]*tau_16_3
		+ REGION4_B[5]*tau_43_3
		+ REGION4_B[6]*tau_110_3;

	return delta * IAPWS97_RHOCRIT;

	/* FIXME iteratively improve vf estimate */
}

double freesteam_region4_rhog_T(double T){

	IAPWS97_APPROXIMATE;

	double tau = 1. - T / IAPWS97_TCRIT;

	double tau_1_6 = pow(tau,1.0/6);

	double tau_2_6 = SQ(tau_1_6);
	double tau_4_6 = SQ(tau_2_6);
	double tau_8_6 = SQ(tau_4_6);
	double tau_16_6 = SQ(tau_8_6);
	double tau_18_6 = tau_16_6 * tau_2_6;
	double tau_37_6 = SQ(tau_18_6) * tau_1_6;
	double tau_71_6 = tau_37_6 * tau_18_6 * tau_16_6;

	double ln_delta =
		  REGION4_C[1]*tau_2_6
		+ REGION4_C[2]*tau_4_6
		+ REGION4_C[3]*tau_8_6
		+ REGION4_C[4]*tau_18_6
		+ REGION4_C[5]*tau_37_6
		+ REGION4_C[6]*tau_71_6;

	return exp(ln_delta) * IAPWS97_RHOCRIT;

	/* FIXME iteratively improve vg estimate */
}


/*------------------------------------------------------------------------------
  INTERPOLATIONS FOR PROPERTIES WITHIN REGION 4
*/

double freesteam_region4_v_Tx(double T, double x){
	double vf, vg;
	if(T < REGION1_TMAX){
		double psat = freesteam_region4_psat_T(T);
		vf = freesteam_region1_v_pT(psat,T);
		vg = freesteam_region2_v_pT(psat,T);
	}else{
		vf = 1./ freesteam_region4_rhof_T(T);
		vg = 1./ freesteam_region4_rhog_T(T);
	}
	return vf + x*(vg - vf);
}

double freesteam_region4_u_Tx(double T, double x){
	double uf, ug;
	if(T < REGION1_TMAX){
		double psat = freesteam_region4_psat_T(T);
		uf = freesteam_region1_u_pT(psat,T);
		ug = freesteam_region2_u_pT(psat,T);
	}else{
		double rhof, rhog;
		rhof = freesteam_region4_rhof_T(T);
		rhog = freesteam_region4_rhog_T(T);
		uf = freesteam_region3_u_rhoT(rhof,T);
		ug = freesteam_region3_u_rhoT(rhog,T);
	}
	return uf + x*(ug - uf);
}

double freesteam_region4_h_Tx(double T, double x){
	double hf, hg;
	if(T < REGION1_TMAX){
		double psat = freesteam_region4_psat_T(T);
		hf = freesteam_region1_h_pT(psat,T);
		hg = freesteam_region2_h_pT(psat,T);
		//fprintf(stderr,"%s: T = %f K, psat = %f MPa, hf = %f kJ/kg, hg = %f kJ/kg\n",__func__,T,psat/1e6,hf/1e3,hg);
	}else{
		double rhof, rhog;
		rhof = freesteam_region4_rhof_T(T);
		rhog = freesteam_region4_rhog_T(T);
		hf = freesteam_region3_h_rhoT(rhof,T);
		hg = freesteam_region3_h_rhoT(rhog,T);
	}
	return hf + x*(hg - hf);
}

double freesteam_region4_s_Tx(double T, double x){
	double sf, sg;
	if(T < REGION1_TMAX){
		double psat = freesteam_region4_psat_T(T);
		sf = freesteam_region1_s_pT(psat,T);
		sg = freesteam_region2_s_pT(psat,T);
	}else{
		double rhof, rhog;
		rhof = freesteam_region4_rhof_T(T);
		rhog = freesteam_region4_rhog_T(T);
		sf = freesteam_region3_s_rhoT(rhof,T);
		sg = freesteam_region3_s_rhoT(rhog,T);
	}
	return sf + x*(sg - sf);
}

double freesteam_region4_cp_Tx(double T, double x){
	double cpf, cpg;
	if(T < REGION1_TMAX){
		double psat = freesteam_region4_psat_T(T);
		cpf = freesteam_region1_cp_pT(psat,T);
		cpg = freesteam_region2_cp_pT(psat,T);
	}else{
		double rhof, rhog;
		rhof = freesteam_region4_rhof_T(T);
		rhog = freesteam_region4_rhog_T(T);
		cpf = freesteam_region3_cp_rhoT(rhof,T);
		cpg = freesteam_region3_cp_rhoT(rhog,T);
	}
	return cpf + x*(cpg - cpf);
}

double freesteam_region4_cv_Tx(double T, double x){
	double cvf, cvg;
	if(T < REGION1_TMAX){
		double psat = freesteam_region4_psat_T(T);
		cvf = freesteam_region1_cv_pT(psat,T);
		cvg = freesteam_region2_cv_pT(psat,T);
	}else{
		double rhof, rhog;
		rhof = freesteam_region4_rhof_T(T);
		rhog = freesteam_region4_rhog_T(T);
		cvf = freesteam_region3_cv_rhoT(rhof,T);
		cvg = freesteam_region3_cv_rhoT(rhog,T);
	}
	return cvf + x*(cvg - cvf);
}

/*------------------------------------------------------------------------------
*/

double freesteam_region4_dpsatdT_T(double T){
	/* calculated this derivative using implicit differentiation of the
	quadratic expression, then derivatives of beta and script-theta */
	double beta = pow(freesteam_region4_psat_T(T)/REGION4_PSTAR, 0.25);
#define N REGION4_N
	double theta = T/REGION4_TSTAR + N[9] / (T/REGION4_TSTAR - N[10]);
	double XBETA = (2.*beta + N[3])*SQ(theta) + (2.*beta*N[1] + N[4])*theta + 2.*N[2]*beta + N[5];
	double XTHETA = (2.*theta + N[1])*SQ(beta) + (2.*N[3]*theta + N[4])*beta + 2.*N[6]*theta + N[7]; 

	double dthetadT = (1 - N[9] / (T/REGION4_TSTAR - N[10]))/REGION4_TSTAR;
	double dbetadtheta = -XTHETA/XBETA;
	double dpdbeta = 4*SQ(beta)*beta*REGION4_PSTAR;
#undef N

	return dpdbeta * dbetadtheta * dthetadT;
}

