/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
