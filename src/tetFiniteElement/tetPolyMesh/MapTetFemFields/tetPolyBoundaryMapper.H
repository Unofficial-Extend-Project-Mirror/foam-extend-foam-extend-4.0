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
    tetPolyBoundaryMapper

Description

SourceFiles
    tetPolyBoundaryMapper.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyBoundaryMapper_H
#define tetPolyBoundaryMapper_H

#include "PtrList.H"
#include "tetPolyPatchMapper.H"
#include "pointMapper.H"
#include "faceMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
              Class tetPolyBoundaryMapper Declaration
\*---------------------------------------------------------------------------*/

class tetPolyBoundaryMapper
:
    public PtrList<tetPolyPatchMapper>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyBoundaryMapper(const tetPolyBoundaryMapper&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyBoundaryMapper&);


public:

    // Constructors

        //- Construct from components
        tetPolyBoundaryMapper
        (
            const tetPolyMesh& mesh,
            const mapPolyMesh& meshMap,
            const pointMapper& pMapper,
            const faceMapper& fMapper
        )
        :
            PtrList<tetPolyPatchMapper>(mesh.boundary().size())
        {
            const tetPolyBoundaryMesh& patches = mesh.boundary();

            forAll (patches, patchI)
            {
                this->set
                (
                    patchI,
                    new tetPolyPatchMapper
                    (
                        patches[patchI],
                        meshMap,
                        pMapper,
                        fMapper
                    )
                );
            }
        }

    // Destructor - default

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
