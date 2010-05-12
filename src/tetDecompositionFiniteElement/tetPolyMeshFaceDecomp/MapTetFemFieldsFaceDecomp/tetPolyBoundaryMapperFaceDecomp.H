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
    tetPolyBoundaryMapperFaceDecomp

Description

SourceFiles
    tetPolyBoundaryMapperFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyBoundaryMapperFaceDecomp_H
#define tetPolyBoundaryMapperFaceDecomp_H

#include "PtrList.H"
#include "tetPolyPatchMapperFaceDecomp.H"
#include "pointMapper.H"
#include "faceMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
              Class tetPolyBoundaryMapperFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyBoundaryMapperFaceDecomp
:
    public PtrList<tetPolyPatchMapperFaceDecomp>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyBoundaryMapperFaceDecomp(const tetPolyBoundaryMapperFaceDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyBoundaryMapperFaceDecomp&);


public:

    // Constructors

        //- Construct from components
        tetPolyBoundaryMapperFaceDecomp
        (
            const tetPolyMeshFaceDecomp& mesh,
            const mapPolyMesh& meshMap,
            const pointMapper& pMapper,
            const faceMapper& fMapper
        )
        :
            PtrList<tetPolyPatchMapperFaceDecomp>(mesh.boundary().size())
        {
            const tetPolyBoundaryMeshFaceDecomp& patches = mesh.boundary();

            forAll (patches, patchI)
            {
                this->set
                (
                    patchI,
                    new tetPolyPatchMapperFaceDecomp
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
