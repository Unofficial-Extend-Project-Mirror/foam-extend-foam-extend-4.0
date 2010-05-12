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
    emptyTetPolyPatchCellDecomp

Description
    Empty plane patch.

SourceFiles
    emptyTetPolyPatchCellDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef emptyTetPolyPatchCellDecomp_H
#define emptyTetPolyPatchCellDecomp_H

#include "faceTetPolyPatchCellDecomp.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class emptyTetPolyPatchCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class emptyTetPolyPatchCellDecomp
:
    public faceTetPolyPatchCellDecomp
{

public:

    //- Runtime type information
    TypeName(emptyPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        emptyTetPolyPatchCellDecomp
        (
            const polyPatch& patch,
            const tetPolyBoundaryMeshCellDecomp& bm
        )
        :
            faceTetPolyPatchCellDecomp(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
