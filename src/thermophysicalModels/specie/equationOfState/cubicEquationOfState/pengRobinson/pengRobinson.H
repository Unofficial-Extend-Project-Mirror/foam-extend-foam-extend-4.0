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
    Foam::pengRobinson

Description
    Peng Robinson equation of state.

    Paper:
    Title: A New Two-Constant Equation of State
    Authors: Ding-Yu Peng and Donald B. Robinson
    Journal: Ind. Eng. Chem., Fundam., Vol. 15, No. 1, 1976


SourceFiles
    pengRobinsonI.H
    pengRobinson.C


Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef pengRobinson_H
#define pengRobinson_H

#include "specie.H"
#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class perfectGas Declaration
\*---------------------------------------------------------------------------*/

class pengRobinson
:
    public specie

{
    // private data
    scalar pcrit_;
    mutable scalar Tcrit_; 
    scalar azentricFactor_;

    //-Peng Robinson factors
    mutable scalar a0_;
    mutable scalar b_;
    mutable scalar n_;

    //CL: pow of constants b_ used in the code e.g. b2_=b*b;
    mutable scalar b2_;
    mutable scalar b3_;
    mutable scalar b4_;
    mutable scalar b5_;
    mutable scalar b6_;

    //CL: rhoMin and rhoMax are only used as boundaries for the bisection methode (see rho function)
    mutable scalar rhoMax_;
    mutable scalar rhoMin_;

    //- Density @STD, initialise after a0, b!
    mutable scalar rhostd_;

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
        inline pengRobinson
        (
	    const specie& sp
	);

        //- Construct from Istream
        pengRobinson(Istream&);

        //- Construct as named copy
        inline pengRobinson(const word& name,const pengRobinson&);

        //- Construct and return a clone
        inline autoPtr<pengRobinson> clone() const;

        // Selector from Istream
        inline static autoPtr<pengRobinson> New(Istream& is);

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

        //Return Peng Robinson factors
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

	inline scalar& b4()const;

	inline scalar& b4();

	inline scalar& b5()const;

	inline scalar& b5();

	inline scalar& b6()const;

	inline scalar& b6();

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
            const  scalar T
        ) const; 

        inline scalar isothermalCompressiblity
        (
            const scalar rho,
            const  scalar T
        ) const; 

        //CL: Used for cv
        inline scalar integral_d2pdT2_dv
        ( 
            const scalar rho,
            const  scalar T
        ) const ; 
	
        //CL: second order derivatives, not Used At The Moment
        inline scalar d2pdv2(const scalar rho,const  scalar T) const;    		
	
        inline scalar d2pdT2(const scalar rho,const  scalar T) const;    		
        
        inline scalar d2pdvdT(const scalar rho,const  scalar T) const;    		
	
        inline scalar d2vdT2(const scalar rho,const  scalar T) const;    
	
        //CL: Used for internal Energy
        inline scalar integral_p_dv(const scalar rho,const  scalar T) const;
        
        //CL: Used for Entropy
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

        inline void operator+=(const pengRobinson&);

        inline void operator*=(const scalar);
        
    // Friend operators

        inline friend pengRobinson operator+
        (
            const pengRobinson&,
            const pengRobinson&
        );

        inline friend pengRobinson operator*
        (
            const scalar s,
            const pengRobinson&
        );

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const pengRobinson&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pengRobinsonI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
