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
    ANY WARRANTY; without even the implied wrranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "oversetMesh.H"
#include "fvcDiv.H"
#include "oversetAdjustPhi.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::oversetMesh::correctNonOrthoFluxes
(
    fvMatrix<Type>& eqn,
    const volVectorField& U
) const
{
    // If the mesh is not orthogonal, get the face flux correction from the
    // matrix, subtract it since it has been already added, make the adjustment
    // and then add it back.
    if
    (
       !mesh().orthogonal()
     && eqn.faceFluxCorrectionPtr()
    )
    {
        typename fvMatrix<Type>::surfaceTypeFieldPtr& nonOrthoFluxPtr =
            eqn.faceFluxCorrectionPtr();

        eqn -= fvc::div(*nonOrthoFluxPtr);
        oversetAdjustPhi(*nonOrthoFluxPtr, U); // Fringe flux adjustment
        eqn += fvc::div(*nonOrthoFluxPtr);

        // Note: no need for globalOversetAdjustPhi since the non-orthogonal
        // correction is switched off for non-coupled boundaries.
    }
}


// ************************************************************************* //
