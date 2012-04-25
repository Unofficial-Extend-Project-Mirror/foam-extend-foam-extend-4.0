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
