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
    Foam::tetherPotential

Description

SourceFiles
    tetherPotential.C
    newTetherPotential.C

\*---------------------------------------------------------------------------*/

#ifndef tetherPotential_H
#define tetherPotential_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class tetherPotential Declaration
\*---------------------------------------------------------------------------*/

class tetherPotential
{

protected:

    // Protected data

        word name_;
        dictionary tetherPotentialProperties_;


    // Private Member Functions

        //- Disallow copy construct
        tetherPotential(const tetherPotential&);

        //- Disallow default bitwise assignment
        void operator=(const tetherPotential&);


public:

    //- Runtime type information
    TypeName("tetherPotential");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            tetherPotential,
            dictionary,
            (
                const word& name,
                const dictionary& tetherPotentialProperties
            ),
            (name, tetherPotentialProperties)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<tetherPotential> New
        (
            const word& name,
            const dictionary& tetherPotentialProperties
        );


    // Constructors

        //- Construct from components
        tetherPotential
        (
            const word& name,
            const dictionary& tetherPotentialProperties
        );


    // Destructor

        virtual ~tetherPotential()
        {}


    // Member Functions

        virtual scalar energy (const vector r) const = 0;

        virtual vector force (const vector r) const = 0;

        const dictionary& tetherPotentialProperties() const
        {
            return tetherPotentialProperties_;
        }

        //- Read tetherPotential dictionary
        virtual bool read(const dictionary& tetherPotentialProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
