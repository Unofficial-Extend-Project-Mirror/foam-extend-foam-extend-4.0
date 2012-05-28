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
    // Static Private data

    //should be read from the fvSolution file where rhoMax and rhoMin values must be defined (for rhoSimpleFoam)
    //HR: Don't know, yet. Let's make these static for starters
    static const scalar rhoMax_;
    static const scalar rhoMin_;

protected:
    // Protected data
    scalar pcrit_;
    scalar Tcrit_; 
    scalar azentricFactor_;

    //-Peng Robinson factors
    scalar n_; 
    scalar a0_;
    mutable scalar b_;


    //CL: Variables to save the values of a, dadT and d2adT2 of the mixture
    //CL: Variables must corrected for changing temperatures
    mutable scalar aSave;
    mutable scalar daSave;
    mutable scalar d2aSave;

    //CL: save the temperature for which the save coefficients (amix,dadTmix,d2adT2mix) are correct
    mutable scalar TSave;

    //- Density @STD, initialise after a0, b!
    mutable scalar rhostd_;


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
	inline scalar rhostd()const;

        //CL: Model coefficient a(T)
	inline scalar a(const  scalar T)const;

        //CL: temperature deriviative of model coefficient a(T)
	inline scalar dadT(const  scalar T)const;

        //CL: second order temperature deriviative of model coefficient a(T)
	inline scalar d2adT2(const  scalar T)const;

	inline scalar a0()const;

	inline scalar b()const;

	inline scalar n()const;

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

        /*
        inline void operator+=(const pengRobinson&);

        inline void operator-=(const pengRobinson&);

        inline void operator*=(const scalar);
        */


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
