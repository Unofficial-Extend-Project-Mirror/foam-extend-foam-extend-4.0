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

Class
    Foam::mixtureSoaveRedlichKwong

Description
    Soave Redlich Kwong equation of state for multiple component mixtures using the van der Waals mixing rule

    Mixing Rule --> see Paper

    Title: Van der eaals mixting rules for cublic equations of state. Applications for supercritical fluids extraction modelling
    Authors: T.Y. Kwak and G.A. Mansoori
    Journal: Chemical Engineering Science, Vol 41, No. 5, pp. 1303-1309, 1986

SourceFiles
    mixtureSoaveRedlichKwongI.H
    mixtureSoaveRedlichKwong.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef mixtureSoaveRedlichKwong_H
#define mixtureSoaveRedlichKwong_H

#include "specie.H"
#include "autoPtr.H"
#include "soaveRedlichKwong.H"
#include "scalarList.H"
#include <vector>
#include "label.H"
#include "List.H"
#include "dictionary.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixtureSoaveRedlichKwong Declaration
\*---------------------------------------------------------------------------*/

class mixtureSoaveRedlichKwong
:
    public soaveRedlichKwong
{

protected:

    //CL: bool used to make sure the model behaves like a single component model when needed
    //CL: this is needed during the construction of the mixture
    mutable bool singleComponent;    
  
    //CL: save the concentrations of each component of the mixture
    //CL: needs to be multiplied by this->W() to get the molar fractions
    mutable std::vector<scalar> weigths;
    
    //CL: stores the values of a(T) of all componentes
    mutable std::vector<scalar> aComponents;

    //CL: stores the values of the derivatives of da/dT of all componentes
    mutable std::vector<scalar> daComponents;

    //CL: stores the values of the second order derivatives of d2a/dT2 of all componentes
    mutable std::vector<scalar> d2aComponents;
    
    //CL: saves a pointer to the pure component classes of the mixture
    mutable std::vector<mixtureSoaveRedlichKwong*> mixtureComponents;
   
    //CL: counts the number of components
    mutable label numOfComp;

    //Protected functions

    //CL: function updates the coefficients (aSave, daSave, d2aSave, b) of the mixture 
    inline void updateModelCoefficients(const  scalar T) const; 

public:

    // Constructors

        //- Construct from components
        //CL: needed for operator*
        inline mixtureSoaveRedlichKwong
        (
            const soaveRedlichKwong& srk,
            label numOfComp,
            std::vector<scalar> weigths,
            std::vector<mixtureSoaveRedlichKwong*> mixtureComponents,
            scalar a0, 
            scalar b,
            scalar Tcrit,
            scalar n,
            scalar rhostd
        );

        //- Construct from components
        //CL: needed for operator+
        inline mixtureSoaveRedlichKwong
        (
            const soaveRedlichKwong& srk,
            label numOfComp,
            std::vector<scalar> weigths,
            std::vector<mixtureSoaveRedlichKwong*> mixtureComponents
        );

        //- Construct from Istream
        mixtureSoaveRedlichKwong(Istream&);

        //- Construct as named copy
        inline mixtureSoaveRedlichKwong(const word& name, const mixtureSoaveRedlichKwong&);

        //- Construct and return a clone
        inline autoPtr<mixtureSoaveRedlichKwong> clone() const;

        // Selector from Istream
        inline static autoPtr<mixtureSoaveRedlichKwong> New(Istream& is);

   // Member functions

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
        ) const;

        //Used for internal Energy
        inline scalar integral_p_dv(const scalar rho, const scalar T) const;  

        // Used for Entropy
        inline scalar integral_dpdT_dv(const scalar rho, const scalar T) const;

        // second order derivatives, not Used At The Moment
        inline scalar d2pdv2(const scalar rho, const scalar T) const;

        inline scalar d2pdT2(const scalar rho, const scalar T) const;

        inline scalar d2pdvdT(const scalar rho, const scalar T) const;

        inline scalar d2vdT2(const scalar rho, const  scalar T) const;

        //- Return density [kg/m^3]
        // rho0 is the starting point of the newton solver used to calculate rho
        inline scalar rho
        (
            const scalar p,
            const scalar T,
            const scalar rho0
        ) const;

        inline scalar rho(const scalar p, const scalar T) const;

        //- Return compressibility drho/dp at T=constant [s^2/m^2]  
        inline scalar psi(const scalar rho, const scalar T) const;

        //- Return compression factor []
        inline scalar Z
        (
            const scalar p,
            const scalar T,
            const scalar rho0
        ) const;


  // Friend operators

        inline friend mixtureSoaveRedlichKwong operator+
        (
            const mixtureSoaveRedlichKwong&,
            const mixtureSoaveRedlichKwong&
        );

        inline friend mixtureSoaveRedlichKwong operator*
        (
            const scalar s,
            const mixtureSoaveRedlichKwong&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const mixtureSoaveRedlichKwong&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mixtureSoaveRedlichKwongI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
