/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    thermalGap

Description
    Thermal properties of a gap.

Author
    Henrik Rusche, Wikki GmbH, All rights reserved.

SourceFiles
    thermalGap.C
    newThermalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef thermalGap_H
#define thermalGap_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalGap Declaration
\*---------------------------------------------------------------------------*/

class thermalGap
{
    // Private data

        //- Name
        const word name_;

        //- Reference to temperature field
        const volScalarField& T_;


    // Private Member Functions

        //- Disallow copy construct
        thermalGap(const thermalGap&);

        //- Disallow default bitwise assignment
        void operator=(const thermalGap&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return T_.mesh();
        }

public:

    //- Runtime type information
    TypeName("thermalGap");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalGap,
            dictionary,
            (
                const word name,
                const volScalarField& T,
                const dictionary& dict
            ),
            (name, T, dict)
        );


    // Selectors

        //- Return a reference to the selected thermal model
        static autoPtr<thermalGap> New
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        thermalGap
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~thermalGap()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return temperture
        const volScalarField& T() const
        {
            return T_;
        };

        //- Modify the resistance
        virtual void modifyResistance(surfaceScalarField&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
