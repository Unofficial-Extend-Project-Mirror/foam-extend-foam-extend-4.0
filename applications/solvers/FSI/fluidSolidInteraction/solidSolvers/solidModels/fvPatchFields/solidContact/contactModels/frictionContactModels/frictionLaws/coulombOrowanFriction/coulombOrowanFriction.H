/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    coulombOrowanFriction

Description
    coulombOrowanFriction friction law
    Essentially the same as Coulomb except the shear traction is limited,
    typically to the material shear yield stress.
    See Computational Contact Mechanics by P.Wriggers

SourceFiles
    coulombOrowanFriction.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef coulombOrowanFriction_H
#define coulombOrowanFriction_H

#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coulombOrowanFriction Declaration
\*---------------------------------------------------------------------------*/

class coulombOrowanFriction
:
    public frictionLaw
{
    // Private data

         //- friction law coefficientdict
        dictionary frictionLawDict_;

        //- Coulomb friction coefficient
        scalar frictionCoeff_;

        //- shear stress limit
        scalar shearStressLimit_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        //coulombOrowanFriction(const coulombOrowanFriction&);

        //- Disallow default bitwise assignment
        void operator=(const coulombOrowanFriction&);

public:

    //- Runtime type information
    TypeName("coulombOrowan");

    // Static data members


    // Constructors

        //- Construct from dictionary
        coulombOrowanFriction
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construct as a copy
        coulombOrowanFriction(const coulombOrowanFriction& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<frictionLaw> clone() const
        {
            return
                autoPtr<frictionLaw>
                (
                    new coulombOrowanFriction(*this)
                );
        }



    // Destructor

        virtual ~coulombOrowanFriction();


    // Member Functions

        //- Return slip traction field
        // the traction that would cause slipping
        // for the given pressure
        virtual scalar slipTraction(const scalar pressure);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
