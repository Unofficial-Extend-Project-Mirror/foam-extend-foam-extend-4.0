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
#include "word.H"
#include "scalar.H"

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


public:
  

    // Constructors

        //- Construct from components
        inline aungierRedlichKwong(
		const specie& sp, 
		scalar pcrit, 
		scalar Tcrit,
		scalar azentricFactor,
	        scalar rhocrit
		);

        //- Construct from Istream
        aungierRedlichKwong(Istream&);

        //- Construct as named copy
        inline aungierRedlichKwong(const word& name, aungierRedlichKwong&);

        //- Construct and return a clone
        inline autoPtr<aungierRedlichKwong> clone() const;

        // Selector from Istream
        inline static autoPtr<aungierRedlichKwong> New(Istream& is);

        //Member Variabels
	    scalar pcrit_;
 	    scalar Tcrit_; 
            scalar rhostd_;
	    scalar azentricFactor_;
 	    scalar rhocrit_;
	    scalar rhoMax_;	//should be read from the fvSolution file where rhoMax and rhoMin values must be define ( for rhoSimpleFoam)
    	    scalar rhoMin_;	

    // Member functions
        inline scalar pcrit() const;
        inline scalar Tcrit() const;
	inline scalar azentricFactor() const;
	inline scalar rhostd()const;
        inline scalar rhocrit() const;
        inline scalar pReturn(const scalar rho, const scalar T) const;

	//-Redlich Kwong factors
     	inline scalar a() const;
	inline scalar b() const;
	inline scalar c() const;
    inline scalar n() const;
        //derivatives 
        inline scalar dpdv(const scalar rho,const  scalar T) const;
        inline scalar dpdT(const scalar rho, const scalar T) const;
        inline scalar dvdT(const scalar rho,const  scalar T) const;
        inline scalar dvdp(const scalar rho, const scalar T) const; 
        inline scalar isobarExpCoef(const scalar rho,const  scalar T) const; 
        inline scalar isothermalCompressiblity(const scalar rho,const  scalar T) const; 
        inline scalar integral_d2pdT2_dv(const scalar rho,const  scalar T) const ; // Used for cv
        inline scalar d2pdv2(const scalar rho,const  scalar T) const;    // not Used At The Moment
	inline scalar d2pdT2(const scalar rho,const  scalar T) const;    // not Used At The Moment
        inline scalar d2pdvdT(const scalar rho,const  scalar T) const;    // not Used At The Moment
	inline scalar d2vdT2(const scalar rho,const  scalar T) const;    // not Used At The Moment
        inline scalar integral_p_dv(const scalar rho,const  scalar T) const;     //Used for internal Energy
        inline scalar integral_dpdT_dv(const scalar rho,const  scalar T) const;  //Used for Entropy

        //- Return density [kg/m^3] // rho0 is the starting point of the newton solver used to calculate rho
        inline scalar rho(const scalar p,const  scalar T,const scalar rho0) const;
        inline scalar rho(const scalar p,const  scalar T) const;

        //- Return compressibility drho/dp [s^2/m^2]   
        inline scalar psi(const scalar rho, const scalar T) const;

        //- Return compression factor []	// rho0 is the starting point of the newton solver used to calculate rho
        inline scalar Z(const scalar p,const  scalar T,const scalar rho0) const;


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
