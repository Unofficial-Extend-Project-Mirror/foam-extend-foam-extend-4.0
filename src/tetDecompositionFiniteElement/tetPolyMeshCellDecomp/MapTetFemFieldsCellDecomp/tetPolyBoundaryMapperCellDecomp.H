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
    tetPolyBoundaryMapperCellDecomp

Description

SourceFiles
    tetPolyBoundaryMapperCellDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyBoundaryMapperCellDecomp_H
#define tetPolyBoundaryMapperCellDecomp_H

#include "PtrList.H"
#include "tetPolyPatchMapperCellDecomp.H"
#include "pointMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
              Class tetPolyBoundaryMapperCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyBoundaryMapperCellDecomp
:
    public PtrList<tetPolyPatchMapperCellDecomp>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyBoundaryMapperCellDecomp(const tetPolyBoundaryMapperCellDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyBoundaryMapperCellDecomp&);


public:

    // Constructors

        //- Construct from components
        tetPolyBoundaryMapperCellDecomp
        (
            const tetPolyMeshCellDecomp& mesh,
            const mapPolyMesh& meshMap
        )
        :
            PtrList<tetPolyPatchMapperCellDecomp>(mesh.boundary().size())
        {
            const tetPolyBoundaryMeshCellDecomp& patches = mesh.boundary();

            forAll (patches, patchI)
            {
                this->set
                (
                    patchI,
                    new tetPolyPatchMapperCellDecomp
                    (
                        patches[patchI],
                        meshMap
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
