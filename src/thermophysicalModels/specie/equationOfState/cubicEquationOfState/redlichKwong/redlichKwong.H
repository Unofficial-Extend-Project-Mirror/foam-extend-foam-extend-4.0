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
    Foam::redlichKwong

Description
    Redlich Kwong equation of state.

SourceFiles
    redlichKwongI.H
    redlichKwong.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef redlichKwong_H
#define redlichKwong_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class redlichKwong Declaration
\*---------------------------------------------------------------------------*/

class redlichKwong
:
    public specie
{


protected:

    // Protected data
    scalar pcrit_;
    scalar Tcrit_; 

    //CL: rhoMin and rhoMax are only used as boundaries for the bisection methode (see rho function)
    scalar rhoMax_;
    scalar rhoMin_;

    //-Redlich Kwong factors
    mutable scalar a_;
    mutable scalar b_;

    //CL: pow of constants b_ used in the code e.g. b2_=b*b;
    mutable scalar b2_;
    mutable scalar b3_;
    mutable scalar b5_;

    //- Density @STD, initialise after a, b!
    mutable scalar rhostd_;

public:
  
    // Constructors

        //- Construct from components
        inline redlichKwong
        (
            const specie& sp
        );

        //- Construct from Istream
        redlichKwong(Istream&);

        //- Construct from dictionary
        //redlichKwong(const dictionary& dict);

        //- Construct as named copy
        inline redlichKwong(const word& name, const redlichKwong&);

        //- Construct and return a clone
        inline autoPtr<redlichKwong> clone() const;

        // Selector from Istream
        inline static autoPtr<redlichKwong> New(Istream& is);

    // Member functions

        inline scalar rhostd() const;

        //Return Redlich Kwong factors
        inline scalar a() const;

        inline scalar b() const;

        inline scalar p(const scalar rho, const scalar T) const;

        //CL: return power of constants b_
	inline scalar b2()const;

	inline scalar b3()const;

	inline scalar b5()const;


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

        inline void operator+=(const redlichKwong&);
/*        
        inline void operator-=(const redlichKwong&);
        inline void operator*=(const scalar);
*/

    // Friend operators

        inline friend redlichKwong operator+
        (
            const redlichKwong&,
            const redlichKwong&
        );

        inline friend redlichKwong operator*
        (
            const scalar s,
            const redlichKwong&
        );

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const redlichKwong&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "redlichKwongI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
