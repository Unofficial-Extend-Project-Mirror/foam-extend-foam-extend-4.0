/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::soaveRedlichKwong

Description
    Soave Redlich Kwong equation of state.

    Paper:
    Title: Equilibrium Constants from a Modified Redlich-Kwong Equation of State
    Authors: G. Soave
    Journal: Chemical Engineering Science, vol. 27(6), 1972, pp. 1197-1203

SourceFiles
    soaveRedlichKwongI.H
    soaveRedlichKwong.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef soaveRedlichKwong_H
#define soaveRedlichKwong_H

#include "specie.H"
#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Soave Readlich Kwong Declaration
\*---------------------------------------------------------------------------*/

class soaveRedlichKwong
:
    public specie

{
    // private data
    scalar pcrit_;
    mutable scalar Tcrit_; 
    scalar azentricFactor_;

    //-Soave Redlich Kwong
    mutable scalar a0_;
    mutable scalar b_;
    mutable scalar n_; 

    //CL: rhoMin and rhoMax are only used as boundaries for the bisection methode (see rho function)
    mutable scalar rhoMin_;
    mutable scalar rhoMax_;

    //- Density @STD, initialise after a0, b!
    mutable scalar rhostd_;

    //CL: pow of constants b_ used in the code e.g. b2_=b*b;
    mutable scalar b2_;
    mutable scalar b3_;
    mutable scalar b5_;

protected:
    //CL: Variables to save the values of a, dadT and d2adT2 of the mixture
    //CL: Variables must corrected for changing temperatures
    mutable scalar aSave;
    mutable scalar daSave;
    mutable scalar d2aSave;

    //CL: save the temperature for which the save coefficients (amix,dadTmix,d2adT2mix) are correct
    mutable scalar TSave;

    //Protected functions
    //CL: function updates the coefficients (aSave, daSave, d2aSave)
    inline void updateModelCoefficients(const  scalar T) const; 

public:
  

    // Constructors

        //- Construct from components
        inline soaveRedlichKwong
        (
            const specie& s
        );

        //- Construct from Istream
        soaveRedlichKwong(Istream&);

        //- Construct as named copy
        inline soaveRedlichKwong(const word& name,const soaveRedlichKwong&);

        //- Construct and return a clone
        inline autoPtr<soaveRedlichKwong> clone() const;

        // Selector from Istream
        inline static autoPtr<soaveRedlichKwong> New(Istream& is);
  
    // Member functions
        inline scalar& rhostd();

        inline scalar& rhostd() const;

        inline scalar& rhoMin();

        inline scalar& rhoMin() const;

        inline scalar& rhoMax();

        inline scalar& rhoMax() const;

        inline scalar& Tcrit();

        inline scalar& Tcrit() const;

        //CL: Model coefficient a(T)
	inline scalar a(const  scalar T)const;

        //CL: temperature deriviative of model coefficient a(T)
	inline scalar dadT(const  scalar T)const;

        //CL: second order temperature deriviative of model coefficient a(T)
	inline scalar d2adT2(const  scalar T)const;

        //Return Soave Redlich Kwong factors
	inline scalar& a0()const;

	inline scalar& a0();

	inline scalar& b()const;

	inline scalar& b();

	inline scalar& n()const;

	inline scalar& n();

        //CL: return power of constants b_
	inline scalar& b2()const;

	inline scalar& b2();

	inline scalar& b3()const;

	inline scalar& b3();

	inline scalar& b5()const;

	inline scalar& b5();

        //CL: Equation of state
        inline scalar p(const scalar rho, const scalar T) const;

        //CL: first order derivatives  
        inline scalar dpdv(const scalar rho,const  scalar T) const;

        inline scalar dpdT(const scalar rho, const scalar T) const;

        inline scalar dvdT(const scalar rho,const  scalar T) const;

        inline scalar dvdp(const scalar rho, const scalar T) const; 

        inline scalar isobarExpCoef
        (
            const scalar rho,
            const scalar T
        ) const; 

        inline scalar isothermalCompressiblity
        (
            const scalar rho,
            const scalar T
        ) const; 

        //CL: Used for cv
        inline scalar integral_d2pdT2_dv(const scalar rho,const  scalar T) const ; 

        //CL: second order derivatives, not Used At The Moment
        inline scalar d2pdv2(const scalar rho,const  scalar T) const;    

	inline scalar d2pdT2(const scalar rho,const  scalar T) const;    

        inline scalar d2pdvdT(const scalar rho,const  scalar T) const;    

	inline scalar d2vdT2(const scalar rho,const  scalar T) const;    

        //CL: Used for internal Energy
        inline scalar integral_p_dv(const scalar rho,const  scalar T) const;    

        //Used for Entropy
        inline scalar integral_dpdT_dv(const scalar rho,const  scalar T) const; 

        //- Return density [kg/m^3] 
        // rho0 is the starting point of the newton solver used to calculate rho
        inline scalar rho
        (
            const scalar p,
            const scalar T,
            const scalar rho0
        ) const;

        inline scalar rho(const scalar p,const  scalar T) const;

        //- Return compressibility drho/dp at T=constant [s^2/m^2] 
        inline scalar psi(const scalar rho, const scalar T) const;

        //- Return compression factor []
        inline scalar Z
        (
            const scalar p,
            const scalar T,
            const scalar rho0
        ) const;

    // Member operators

        inline void operator+=(const soaveRedlichKwong&);

        inline void operator*=(const scalar);

    // Friend operators

        inline friend soaveRedlichKwong operator+
        (
            const soaveRedlichKwong&,
            const soaveRedlichKwong&
        );

        inline friend soaveRedlichKwong operator*
        (
            const scalar s,
            const soaveRedlichKwong&
        );

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const soaveRedlichKwong&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "soaveRedlichKwongI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
