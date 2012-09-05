/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    IAPWS-IF97 (water) based thermodynamic class. Water properties calculated by freeSteam.    
    
    This code connects OpenFoam with freeSteam and provides the basic functions needed in OpenFOAM

    For more information about freeSteam and its authors have a look @ http://freesteam.sourceforge.net/example.php

SourceFiles
  IAPWS-IF97.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany
\*---------------------------------------------------------------------------*/
#ifndef IAPWSIF97_H
#define IAPWSIF97_H

#include "basicPsiThermo.H"
#include "steam.H"

#ifdef __cplusplus
	#define EXTERN extern "C"
#else
	#define EXTERN extern
#endif 


//CL: 
EXTERN double freesteam_p(SteamState S);
EXTERN double freesteam_T(SteamState S);
EXTERN double freesteam_rho(SteamState S);
EXTERN double freesteam_v(SteamState S);
EXTERN double freesteam_u(SteamState S);
EXTERN double freesteam_h(SteamState S);
EXTERN double freesteam_s(SteamState S);
EXTERN double freesteam_cp(SteamState S);
EXTERN double freesteam_cv(SteamState S);
EXTERN double freesteam_w(SteamState S);
EXTERN double freesteam_x(SteamState S);
EXTERN double freesteam_mu(SteamState S);
EXTERN double freesteam_k(SteamState S);

//CL: getting SteamState for two given properties e.g. pressure and temperatur
EXTERN SteamState freesteam_set_pv(double,double);
EXTERN SteamState freesteam_set_pu(double,double);
EXTERN SteamState freesteam_set_pT(double,double);
EXTERN SteamState freesteam_set_ph(double,double);

//CL: getting region of the SteamState
EXTERN int freesteam_region(SteamState);

//CL: transport properties
EXTERN double freesteam_mu_rhoT(double,double);
EXTERN double freesteam_k_rhoT(double,double);

//CL: Region 1 --> see region1.h (freesteam)
EXTERN double freesteam_region1_v_pT(double,double);
EXTERN double freesteam_region1_h_pT(double,double);
EXTERN double freesteam_region1_kappaT_pT(double,double);
EXTERN double freesteam_region1_alphav_pT(double,double);
EXTERN double freesteam_region1_cp_pT(double,double);
EXTERN double freesteam_region1_u_pT(double,double);
EXTERN double freesteam_region1_s_pT(double,double);
EXTERN double freesteam_region1_cv_pT(double,double);

//CL: Region 2 --> see region2.h (freesteam)
EXTERN double freesteam_region2_v_pT(double,double);
EXTERN double freesteam_region2_u_pT(double,double);
EXTERN double freesteam_region2_s_pT(double,double);
EXTERN double freesteam_region2_h_pT(double,double);
EXTERN double freesteam_region2_cp_pT(double,double);
EXTERN double freesteam_region2_cv_pT(double,double);
EXTERN double freesteam_region2_alphav_pT(double,double);
EXTERN double freesteam_region2_kappaT_pT(double,double);

//CL: Region 3 --> see region3.h (freesteam)
EXTERN double freesteam_region3_p_rhoT(double,double);
EXTERN double freesteam_region3_u_rhoT(double,double);
EXTERN double freesteam_region3_s_rhoT(double,double);
EXTERN double freesteam_region3_h_rhoT(double,double);
EXTERN double freesteam_region3_cp_rhoT(double,double);
EXTERN double freesteam_region3_cv_rhoT(double,double);
EXTERN double freesteam_region3_alphap_rhoT(double,double);
EXTERN double freesteam_region3_betap_rhoT(double,double);

//CL: Region 4 --> see region4.h (freesteam)
EXTERN double freesteam_region4_psat_T(double);
EXTERN double freesteam_region4_Tsat_p(double);
EXTERN double freesteam_region4_rhof_T(double);
EXTERN double freesteam_region4_rhog_T(double);
EXTERN double freesteam_region4_v_Tx(double,double);
EXTERN double freesteam_region4_u_Tx(double,double);
EXTERN double freesteam_region4_h_Tx(double,double);
EXTERN double freesteam_region4_s_Tx(double,double);
EXTERN double freesteam_region4_cp_Tx(double,double);
EXTERN double freesteam_region4_cv_Tx(double,double);
EXTERN double freesteam_region4_dpsatdT_T(double);

namespace Foam
{
    //CL: Functions to caluculate all fluid properties
    void calculateProperties_h
    (
        SteamState S,
        scalar &rho,
        scalar &h,
        scalar &T,
        scalar &p,
        scalar &psi,
        scalar &drhodh,
        scalar &mu,
        scalar &alpha,
        scalar &x
    );

    //CL: This functions returns all (minimal) needed propeties (p,T,h,rho,psi,drhodh,mu and alpha) for given p and T
    void calculateProperties_pT
    (
        scalar &p,
        scalar &T, 
        scalar &h, 
        scalar &rho, 
        scalar &psi, 
        scalar &drhodh, 
        scalar &mu, 
        scalar &alpha
    );

    //CL: This function returns the same values as the function above for given p and T
    //CL: Additionally, the vapor mass fraction x is return
    //CL: NOTE: This function is only included to have the possibility to update x at the fixedValue (Temperature) BC
    //CL: can only return x=0 and x=1 because it is not possible to describe the vapour dome with p and T 
    void calculateProperties_pT
    (
        scalar &p, 
        scalar &T, 
        scalar &h, 
        scalar &rho, 
        scalar &psi, 
        scalar &drhodh, 
        scalar &mu, 
        scalar &alpha, 
        scalar &x
    );

    //CL: This functions returns all (minimal) needed properties (p,T,h,rho,psi,drhodh,mu and alpha) for given p and h
    void calculateProperties_ph
    (
        scalar &p, 
        scalar &h, 
        scalar &T, 
        scalar &rho, 
        scalar &psi, 
        scalar &drhodh, 
        scalar &mu, 
        scalar &alpha
    );

    //CL: This function returns the same values as the function above for given p and h
    //CL: Additionally, the vapor mass fraction x is return
    void calculateProperties_ph
    (
        scalar &p, 
        scalar &h, 
        scalar &T, 
        scalar &rho, 
        scalar &psi, 
        scalar &drhodh, 
        scalar &mu, 
        scalar &alpha, 
        scalar &x
    );


    //CL: Return density for given pT or ph;
    scalar rho_pT(scalar p,scalar T);
    scalar rho_ph(scalar p,scalar h);

    //CL: Return cp for given pT or ph;
    scalar cp_pT(scalar p,scalar T);
    scalar cp_ph(scalar p,scalar h);

    //CL: Return cv for given pT or ph;
    scalar cv_pT(scalar p,scalar T);
    scalar cv_ph(scalar p,scalar h);

    //CL: Return enthalpy for given pT;
    scalar h_pT(scalar p,scalar T);

    //CL: Return temperature for given ph;
    scalar T_ph(scalar p,scalar T);

    //CL: Return psiH=(drho/dp)_h=constant for given pT or ph;
    scalar psiH_pT(scalar p,scalar T);
    scalar psiH_ph(scalar p,scalar h);
    scalar psiH(SteamState S);

    //CL: Return drhodh=(drho/dh)_p=constant for given pT or ph;
    scalar drhodh_pT(scalar p,scalar T);
    scalar drhodh_ph(scalar p,scalar h);
    scalar drhodh(SteamState S);

}

#endif //IAPWSIF97_C_
