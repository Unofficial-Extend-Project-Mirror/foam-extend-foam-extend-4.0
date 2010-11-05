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
    fluxCorrector

Description
    Virtual base class that deals with flux-correction after topo-changes.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    fluxCorrector.C

\*---------------------------------------------------------------------------*/

#ifndef fluxCorrector_H
#define fluxCorrector_H

#include "fvMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations

/*---------------------------------------------------------------------------*\
                         Class fluxCorrector Declaration
\*---------------------------------------------------------------------------*/

class fluxCorrector
{
    // Private data

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Reference to dictionary
        const dictionary& dict_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fluxCorrector(const fluxCorrector&);

        //- Disallow default bitwise assignment
        void operator=(const fluxCorrector&);

public:

    //- Runtime type information
    TypeName("fluxCorrector");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxCorrector,
            mesh,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Selectors

        //- Select constructed from fvMesh
        static autoPtr<fluxCorrector> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from fvMesh and dictionary
        fluxCorrector(const fvMesh& mesh, const dictionary& dict)
        :
            mesh_(mesh),
            dict_(dict)
        {}


    // Destructor

        virtual ~fluxCorrector()
        {}


    // Member Functions

        //- Return reference to mesh
        const fvMesh& mesh() const;

        //- Return reference to dictionary
        const dictionary& dict() const;

        //- Is flux-correction required?
        virtual bool required() const;

        //- Interpolate fluxes to a specified list of faces
        virtual void interpolateFluxes(const labelList& faces) const;

        //- Update fluxes in the registry
        virtual void updateFluxes() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "fluxCorrector.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
