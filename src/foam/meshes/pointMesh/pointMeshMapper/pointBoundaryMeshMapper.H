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
    Foam::pointBoundaryMeshMapper

Description
    Foam::pointBoundaryMeshMapper

\*---------------------------------------------------------------------------*/

#ifndef pointBoundaryMeshMapper_H
#define pointBoundaryMeshMapper_H

#include "PtrList.H"
#include "pointPatchMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class pointBoundaryMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class pointBoundaryMeshMapper
:
    public PtrList<pointPatchMapper>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointBoundaryMeshMapper(const pointBoundaryMeshMapper&);

        //- Disallow default bitwise assignment
        void operator=(const pointBoundaryMeshMapper&);


public:

    // Constructors

        //- Construct from components
        pointBoundaryMeshMapper
        (
            const pointMesh& mesh,
            const pointMapper& pointMap,
            const mapPolyMesh& mpm
        )
        :
            PtrList<pointPatchMapper>(mesh.boundary().size())
        {
            const pointBoundaryMesh& patches = mesh.boundary();

            forAll (patches, patchI)
            {
                set
                (
                    patchI,
                    new pointPatchMapper
                    (
                        patches[patchI],
                        pointMap,
                        mpm
                    )
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
