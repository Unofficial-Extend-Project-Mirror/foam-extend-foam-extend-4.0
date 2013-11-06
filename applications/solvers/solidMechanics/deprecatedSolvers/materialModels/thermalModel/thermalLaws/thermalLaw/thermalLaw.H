/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008 Hrvoje Jasak
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
    thermalLaw

Description
    Thermal material properties for solids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    thermalLaw.C
    newThermalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef thermalLaw_H
#define thermalLaw_H

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
                         Class thermalLaw Declaration
\*---------------------------------------------------------------------------*/

class thermalLaw
{
    // Private data

        //- Name
        const word name_;

        //- Reference to stress field
        const volScalarField& T_;


    // Private Member Functions

        //- Disallow copy construct
        thermalLaw(const thermalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const thermalLaw&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return T_.mesh();
        }

public:

    //- Runtime type information
    TypeName("thermalLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalLaw,
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
        static autoPtr<thermalLaw> New
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        thermalLaw
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~thermalLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }


        //- Return specific heat capacity
        virtual tmp<volScalarField> C() const = 0;

        //- Return thermal conductivity
        virtual tmp<volScalarField> k() const = 0;

        //- Return thermal expansion coefficient
        virtual tmp<volScalarField> alpha() const = 0;

        //- Return reference temperature
        virtual tmp<volScalarField> T0() const = 0;

        //- Correct the rheological model
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
