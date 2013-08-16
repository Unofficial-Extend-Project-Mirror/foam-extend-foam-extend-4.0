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
    Foam::mixturePengRobinson

Description
    Peng Robinson equation of state for multiple component mixtures using the van der Waals mixing rule

    Mixing Rule --> see Paper

    Title: Van der waals mixing rules for cublic equations of state. Applications for supercritical fluids extraction modelling
    Authors: T.Y. Kwak and G.A. Mansoori
    Journal: Chemical Engineering Science, Vol 41, No. 5, pp. 1303-1309, 1986

SourceFiles
    mixturePengRobinsonI.H
    mixturePengRobinson.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef mixturePengRobinson_H
#define mixturePengRobinson_H

#include "specie.H"
#include "autoPtr.H"
#include "pengRobinson.H"
#include "scalarList.H"
#include "DynamicList.H"
#include "label.H"
#include "List.H"
#include "dictionary.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixturePengRobinson Declaration
\*---------------------------------------------------------------------------*/

class mixturePengRobinson
:
    public pengRobinson
{

protected:

    //CL: bool used to make sure the model behaves like a single component model when needed
    //CL: this is needed during the construction of the mixture
    mutable bool singleComponent;    
  
    //CL: save the concentrations of each component of the mixture
    //CL: needs to be multiplied by this->W() to get the molar fractions
    mutable DynamicList<scalar> weigths;
    
    //CL: saves a pointer to the pure component classes of the mixture
    mutable DynamicList<mixturePengRobinson*> mixtureComponents;
   
    //CL: counts the number of components
    mutable label numOfComp;

    //Protected functions

    //CL: function updates the coefficients (aSave, daSave, d2aSave, b) of the mixture 
    //CL: this is the function with the mixing rule
    inline void updateModelCoefficients(const  scalar T) const; 

    //CL: Variables used in real gas mixture correction

        //CL: If true, the real gas mixture correction is used
        mutable bool  realMixtureCorr_;
 
        //CL: number of mixture components, needed to calculate the mixture correction factors needed
        //CL: do not mistake this variable with numOfComp,
        //CL: numOfComp is a counter to counts the number of components while the mixture is constructed
        mutable label nCom_;
 
        //CL: stores real mixture correction coefficients
        mutable DynamicList<scalar> realMixtureCorrCoef_;


public:

    // Constructors

        //- Construct from components
        //CL: needed for operator*
        inline mixturePengRobinson
        (
            const pengRobinson& pr,
            label numOfComp,
            DynamicList<scalar> weigths,
            DynamicList<mixturePengRobinson*> mixtureComponents,
            scalar a0, 
            scalar b,
            scalar Tcrit,
            scalar n,
            scalar rhostd,
            scalar rhoMin,
            scalar rhoMax
        );

        //- Construct from components
        //CL: needed for operator+
        inline mixturePengRobinson
        (
            const pengRobinson& pr,
            label numOfComp,
            DynamicList<scalar> weigths,
            DynamicList<mixturePengRobinson*> mixtureComponents,
            scalar rhoMin,
            scalar rhoMax
        );

        //- Construct from Istream
        mixturePengRobinson(Istream&);

        //- Construct from dictionary
        //mixturePengRobinson(const dictionary& dict);

        //- Construct as named copy
        inline mixturePengRobinson(const word& name, const mixturePengRobinson&);

        //- Construct and return a clone
        inline autoPtr<mixturePengRobinson> clone() const;

        // Selector from Istream
        inline static autoPtr<mixturePengRobinson> New(Istream& is);

   // Member functions

        //CL: mixture coefficience
        inline scalar amix(const scalar T, const label i, const  label j) const;
                 
        inline scalar dadTmix(const scalar T, const label i, const  label j) const;

        inline scalar d2adT2mix(const scalar T, const label i, const  label j) const;

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

        // I-O

            //- Write to Ostream
            //void write(Ostream& os) const;

  // Member operators

        inline void operator+=(const mixturePengRobinson&);

  // Friend operators

        inline friend mixturePengRobinson operator+
        (
            const mixturePengRobinson&,
            const mixturePengRobinson&
        );

        inline friend mixturePengRobinson operator*
        (
            const scalar s,
            const mixturePengRobinson&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const mixturePengRobinson&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mixturePengRobinsonI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
