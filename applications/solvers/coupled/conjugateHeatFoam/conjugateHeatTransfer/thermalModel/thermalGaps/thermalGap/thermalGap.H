/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 H. Rusche
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
