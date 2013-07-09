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
    thermalSource

Description
    A thermal source.

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved.

SourceFiles
    thermalSource.C
    newThermalSource.C

\*---------------------------------------------------------------------------*/

#ifndef thermalSource_H
#define thermalSource_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalSource Declaration
\*---------------------------------------------------------------------------*/

class thermalSource
{
    // Private data

        //- Name
        const word name_;

        //- Reference to temperature field
        const volScalarField& T_;


    // Private Member Functions

        //- Disallow copy construct
        thermalSource(const thermalSource&);

        //- Disallow default bitwise assignment
        void operator=(const thermalSource&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return T_.mesh();
        }

public:

    //- Runtime type information
    TypeName("thermalSource");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalSource,
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
        static autoPtr<thermalSource> New
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        thermalSource
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~thermalSource()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Add source term
        virtual void addSource(volScalarField& source) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
