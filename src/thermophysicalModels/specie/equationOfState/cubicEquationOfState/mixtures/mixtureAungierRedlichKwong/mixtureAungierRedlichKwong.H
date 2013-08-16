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
    Foam::mixtureAungierRedlichKwong

Description
    Mixture Aungier Redlich Kwong equation of state.

    Mixture based on a pseudo critical point approach. 
    
    For further information, see:
    BOOK Title: The Properties of Gases And Liquids, Fifth Edition, McGraw-Hill, Chapter 5
    Authors: B. Poling; J. Prausnitz; J. O'Connell


SourceFiles
    mixtureAungierRedlichKwongI.H
    mixtureAungierRedlichKwong.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef mixtureAungierRedlichKwong_H
#define mixtureAungierRedlichKwong_H

#include "specie.H"
#include "aungierRedlichKwong.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixtureAungierRedlichKwong Declaration
\*---------------------------------------------------------------------------*/

class mixtureAungierRedlichKwong
:
    public aungierRedlichKwong

{

protected:

    // data at the critical point
    scalar Zcrit_;
        
    //CL. molar volume at the critical point
    scalar Vcrit_;

    //CL: save the concentrations of each component of the mixture
    //CL: needs to be multiplied by this->W() to get the molar fractions
    //mutable DynamicList<scalar> weigths;
    
    //CL: saves a pointer to the pure component classes of the mixture
   // mutable DynamicList<mixtureAungierRedlichKwong*> mixtureComponents;


public:
  

    // Constructors

        //- Construct from components
        inline mixtureAungierRedlichKwong
        (
            const aungierRedlichKwong& sp, 
            scalar Tcrit,
            scalar azentricFactor,
            scalar Vcrit,
            scalar Zcrit,
            scalar rhoMin,
            scalar rhoMax
        );

        //- Construct from Istream
        mixtureAungierRedlichKwong(Istream&);

        //- Construct from dictionary
        //mixtureAungierRedlichKwong(const dictionary& dict);

        //- Construct as named copy
        inline mixtureAungierRedlichKwong(const word& name, const mixtureAungierRedlichKwong&);

        //- Construct and return a clone
        inline autoPtr<mixtureAungierRedlichKwong> clone() const;

        // Selector from Istream
        inline static autoPtr<mixtureAungierRedlichKwong> New(Istream& is);

   // Member functions

        // I-O

            //- Write to Ostream
            //void write(Ostream& os) const;


  // Member operators

        inline void operator+=(const mixtureAungierRedlichKwong&);

  // Friend operators

        inline friend mixtureAungierRedlichKwong operator+
        (
            const mixtureAungierRedlichKwong&,
            const mixtureAungierRedlichKwong&
        );

        inline friend mixtureAungierRedlichKwong operator*
        (
            const scalar s,
            const mixtureAungierRedlichKwong&
        );

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const mixtureAungierRedlichKwong&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mixtureAungierRedlichKwongI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
