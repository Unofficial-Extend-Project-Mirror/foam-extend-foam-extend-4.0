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

/*
 Based on IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance
 */

#define FREESTEAM_BUILDING_LIB
#include "viscosity.h"

static double mu0(double tau);
static double mu1(double del, double tau);

#define VISCOSITY_MUSTAR 1.0e-6 /* Pa-s */

#include <math.h>

static double mu0(double tau){
    // viscosity in the dilute-gas limit
    const double H[4] = {1.67752, 2.20462, 0.6366564, -0.241605};
    int i;
    double sum = 0;
    for (i = 0; i < 4; i++){
        sum += H[i] * ipow(tau, i) ;
    }
    return 100.0 / (sqrt(tau) * sum);
}



static double mu1(double del, double tau){
    // contribution to viscosity due to finite density
    const double H[6][7] = {
        { 5.20094E-1, 2.22531E-1, -2.81378E-1, 1.61913E-1, -3.25372E-2, 0.0,         0.0},
        { 8.50895E-2, 9.99115E-1, -9.06851E-1, 2.57399E-1,  0.0,        0.0,         0.0},
        {-1.08374,    1.88797,    -7.72479E-1, 0.0,         0.0,        0.0,         0.0},
        {-2.89555E-1, 1.26613,    -4.89837E-1, 0.0,         6.98452E-2, 0.0,        -4.35673E-3},
        { 0.0,        0.0,        -2.57040E-1, 0.0,         0.0,        8.72102E-3,  0.0},
        { 0.0,        1.20573E-1,  0.0,        0.0,         0.0,        0.0,        -5.93264E-4}
    };
	
    int i, j;
    double sum = 0;
	double tau1 = 0;
    for (i = 0; i < 6; i++){
		tau1 = ipow(tau - 1, i);
        for (j = 0; j < 7; j++){
			if(0==H[i][j])continue;
            sum += H[i][j] * tau1 * ipow(del - 1, j);
        }
    }
    return exp(del * sum);
}

double freesteam_mu_rhoT(double rho, double T){
	double del = rho / IAPWS97_RHOCRIT;
	double tau = IAPWS97_TCRIT / T;

    const int mu2 = 1; // critical enhancement to viscosity not implemented for IF-97, set to 1
	return VISCOSITY_MUSTAR * mu0(tau) * mu1(del,tau) * mu2;
}

