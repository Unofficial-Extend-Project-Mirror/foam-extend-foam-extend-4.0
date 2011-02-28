/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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

Class
    Foam::aungierRedlichKwong

Description
    Aungier Redlich Kwong equation of state.

SourceFiles
    aungierRedlichKwongI.H
    aungierRedlichKwong.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef aungierRedlichKwong_H
#define aungierRedlichKwong_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class aungierRedlichKwong Declaration
\*---------------------------------------------------------------------------*/

class aungierRedlichKwong
:
    public specie

{
    // Private data
        scalar pcrit_;
        scalar Tcrit_; 
        scalar azentricFactor_;
        scalar rhocrit_;

    //Aungier Redlich Kwong factors
        scalar a_;
        scalar b_;
        scalar c_;
        scalar n_;

    //Density @STD, initialise after a, b!
        scalar rhostd_;

    // Static Private data

    //should be read from the fvSolution file where rhoMax and rhoMin values must be define ( for rhoSimpleFoam)
        static const scalar rhoMax_;
        static const scalar rhoMin_; 


public:
  

    // Constructors

        //- Construct from components
        inline aungierRedlichKwong
        (
            const specie& sp, 
            scalar pcrit, 
            scalar Tcrit,
            scalar azentricFactor,
            scalar rhocrit
        );

        //- Construct from Istream
        aungierRedlichKwong(Istream&);

        //- Construct as named copy
        inline aungierRedlichKwong(const word& name, const aungierRedlichKwong&);

        //- Construct and return a clone
        inline autoPtr<aungierRedlichKwong> clone() const;

        // Selector from Istream
        inline static autoPtr<aungierRedlichKwong> New(Istream& is);

   // Member functions

	inline scalar rhostd() const;

        inline scalar p(const scalar rho, const scalar T) const;

        //first order derivatives 
        inline scalar dpdv(const scalar rho, const  scalar T) const;

        inline scalar dpdT(const scalar rho, const scalar T) const;

        inline scalar dvdT(const scalar rho, const  scalar T) const;

        inline scalar dvdp(const scalar rho, const scalar T) const; 

        inline scalar isobarExpCoef(const scalar rho, const  scalar T) const; 

        inline scalar isothermalCompressiblity
        (
            const scalar rho,
            const scalar T
        ) const; 
        
        // Used for cv
        inline scalar integral_d2pdT2_dv
        (
            const scalar rho,
            const scalar T
        ) const ; 

        // second order derivatives, not Used At The Moment
        inline scalar d2pdv2(const scalar rho, const  scalar T) const;    

	inline scalar d2pdT2(const scalar rho,   scalar T) const;    

        inline scalar d2pdvdT(const scalar rho, const  scalar T) const;    

	inline scalar d2vdT2(const scalar rho, const  scalar T) const;    

        //Used for internal Energy
        inline scalar integral_p_dv(const scalar rho, const  scalar T) const;

        //Used for Entropy
        inline scalar integral_dpdT_dv(const scalar rho, const  scalar T) const;  

        //- Return density [kg/m^3] 
        // rho0 is the starting point of the newton solver used to calculate rho
        inline scalar rho
	(
            const scalar p,
            const scalar T,
            const scalar rho0
        ) const;
 
        inline scalar rho(const scalar p, const  scalar T) const;

        //- Return compressibility drho/dp at T=constant [s^2/m^2]
        inline scalar psi(const scalar rho, const scalar T) const;

        //- Return compression factor  []
        inline scalar Z
        (
            const scalar p,
            const  scalar T, 
            const scalar rho0
        ) const;


    // Member operators

        inline void operator+=(const aungierRedlichKwong&);
        inline void operator-=(const aungierRedlichKwong&);

        inline void operator*=(const scalar);


 /*   // Friend operators

        inline friend aungierRedlichKwong operator+
        (
            const aungierRedlichKwong&,
            const aungierRedlichKwong&
        );

        inline friend aungierRedlichKwong operator-
        (
            const aungierRedlichKwong&,
            const aungierRedlichKwong&
        );

        inline friend aungierRedlichKwong operator*
        (
            const scalar s,
            const aungierRedlichKwong&
        );

        inline friend aungierRedlichKwong operator==
        (
            const aungierRedlichKwong&,
            const aungierRedlichKwong&
        );
*/

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const aungierRedlichKwong&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "aungierRedlichKwongI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
