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
    Foam::faBoundaryMeshMapper

Description
    Foam::faBoundaryMeshMapper

\*---------------------------------------------------------------------------*/

#ifndef faBoundaryMeshMapper_H
#define faBoundaryMeshMapper_H

#include "PtrList.H"
#include "faPatchMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class faBoundaryMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class faBoundaryMeshMapper
:
    public PtrList<faPatchMapper>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        faBoundaryMeshMapper(const faBoundaryMeshMapper&);

        //- Disallow default bitwise assignment
        void operator=(const faBoundaryMeshMapper&);


public:

    // Constructors

        //- Construct from components
        faBoundaryMeshMapper
        (
            const faMesh& mesh,
            const mapPolyMesh& mpm
        )
        :
            PtrList<faPatchMapper>(mesh.boundary().size())
        {
            const faBoundaryMesh& patches = mesh.boundary();

            forAll (patches, patchI)
            {
                set
                (
                    patchI,
                    new faPatchMapper
                    (
                        patches[patchI],
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
