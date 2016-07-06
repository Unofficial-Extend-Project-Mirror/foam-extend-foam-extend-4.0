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

Class
    Foam::aungierRedlichKwong

Description
    Aungier Redlich Kwong equation of state.

    Paper:
    Title: A Fast, Accurate Real Gas Equation of State for Fluid Dynamic Analysis Applications
    Authors: R. H. Aungier
    Journal: Journal of Fluids Engineering, Vol.117, 1995

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
    // private data
    //CL: data at critical point
    scalar pcrit_;
    scalar Tcrit_;
    scalar rhocrit_;
    scalar azentricFactor_;

    //Aungier Redlich Kwong factors
    scalar a0_;
    scalar b_;
    scalar c_;
    scalar n_;

    //CL: pow of constants (b_, c_) used in the code e.g. b2_=b*b;
    scalar b2_;
    scalar b3_;
    scalar b4_;
    scalar b5_;
    scalar c2_;

    //CL: rhoMin and rhoMax are only used as boundaries for the bisection methode (see rho function)
    scalar rhoMin_;
    scalar rhoMax_;

    //CL: Variables to save the values of a, dadT and d2adT2 of the mixture
    //CL: Variables must corrected for changing temperatures
    mutable scalar aSave;
    mutable scalar daSave;
    mutable scalar d2aSave;

    //CL: save the temperature for which the save coefficients (aSave, daSave, d2aSave) are correct
    mutable scalar TSave;

    //- Density @STD, initialise after a0, b, c, n!
    scalar rhostd_;

    //Protected functions
    //CL: function updates the coefficients (aSave, daSave, d2aSave)
    inline void updateModelCoefficients(const  scalar T) const;

public:

    // Constructors

        //- Construct from components
        inline aungierRedlichKwong
        (
            const specie& sp
        );

        //- Construct from Istream
        aungierRedlichKwong(Istream&);

        //- Construct from dictionary
        //aungierRedlichKwong(const dictionary& dict);

        //- Construct as named copy
        inline aungierRedlichKwong(const word& name, const aungierRedlichKwong&);

        //- Construct and return a clone
        inline autoPtr<aungierRedlichKwong> clone() const;

        // Selector from Istream
        inline static autoPtr<aungierRedlichKwong> New(Istream& is);

   // Member functions

        //Return Aungier Redlich Kwong factors
        inline scalar a0() const;

        inline scalar b() const;

        inline scalar c() const;

        inline scalar n() const;

        inline scalar rhostd() const;

        inline scalar rhoMin() const;

        inline scalar rhoMax() const;

        inline scalar Tcrit() const;

        inline scalar pcrit() const;

        inline scalar rhocrit() const;

        inline scalar azentricFactor() const;

        //CL: Model coefficient a(T)
        inline scalar a(const scalar T) const;

        //CL: temperature derivative of model coefficient a(T)
        inline scalar dadT(const scalar T) const;

        //CL: second order temperature derivative of model coefficient a(T)
        inline scalar d2adT2(const scalar T) const;

        //CL: Equation of state
        inline scalar p(const scalar rho, const scalar T) const;

        //CL: first order derivatives
        inline scalar dpdv(const scalar rho, const scalar T) const;

        inline scalar dpdT(const scalar rho, const scalar T) const;

        inline scalar dvdT(const scalar rho, const scalar T) const;

        inline scalar dvdp(const scalar rho, const scalar T) const;

        inline scalar isobarExpCoef(const scalar rho, const scalar T) const;

        inline scalar isothermalCompressiblity
        (
            const scalar rho,
            const scalar T
        ) const;

        //CL: Used for cv
        inline scalar integral_d2pdT2_dv
        (
            const scalar rho,
            const scalar T
        ) const;

        //CL: second order derivatives, not Used At The Moment
        inline scalar d2pdv2(const scalar rho, const scalar T) const;

        inline scalar d2pdT2(const scalar rho, const scalar T) const;

        inline scalar d2pdvdT(const scalar rho, const scalar T) const;

        inline scalar d2vdT2(const scalar rho, const scalar T) const;

        //CL: Used for internal Energy
        inline scalar integral_p_dv
        (
            const scalar rho,
            const scalar T
        ) const;

        //CL: Used for Entropy
        inline scalar integral_dpdT_dv
        (
            const scalar rho,
            const scalar T
        ) const;

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

        inline void operator+=(const aungierRedlichKwong&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend aungierRedlichKwong operator+
        (
            const aungierRedlichKwong&,
            const aungierRedlichKwong&
        );

        inline friend aungierRedlichKwong operator*
        (
            const scalar s,
            const aungierRedlichKwong&
        );


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
