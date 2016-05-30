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


Description:

Includes the definition of a struct from freestream
Code copied from freestream --> steam.h

\*---------------------------------------------------------------------------*/

#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN extern
#endif

typedef struct SteamState_R1_struct{
    double p, T;
} SteamState_R1;

typedef struct SteamState_R2_struct{
    double p, T;
} SteamState_R2;

typedef struct SteamState_R3_struct{
    double rho, T;
} SteamState_R3;

typedef struct SteamState_R4_struct{
    double T, x;
} SteamState_R4;

typedef struct SteamState_struct{
    char region;
    union{
        SteamState_R1 R1;
        SteamState_R2 R2;
        SteamState_R3 R3;
        SteamState_R4 R4;
    };
} SteamState;

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

