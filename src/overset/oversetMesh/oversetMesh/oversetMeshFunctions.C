/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
