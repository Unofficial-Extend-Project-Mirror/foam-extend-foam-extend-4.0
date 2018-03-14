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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

Description
    Specilalisation of manipulateMatrix member function for scalars needed to
    ensure conservation across GGI patches that are partially covered if the
    bridge overlap switch is on.

\*---------------------------------------------------------------------------*/

#include "ggiFvPatchField.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<>
void Foam::ggiFvPatchField<scalar>::manipulateMatrix(fvMatrix<scalar>& matrix)
{
    Info<< "IN MANIPULATE MATRIX!" << endl;

    // Conservative treatment for bridged overlap
    if (ggiPatch_.bridgeOverlap())
    {
        // Get this patch and shadow patch index
        const label patchI = this->patch().index();
        const label sPatchI = ggiPatch_.shadowIndex();

        // Get all matrix coefficients
        scalarField& thisIC = matrix.internalCoeffs()[patchI];
        scalarField& thisBC = matrix.boundaryCoeffs()[patchI];
        scalarField& shadowIC = matrix.internalCoeffs()[sPatchI];
        scalarField& shadowBC = matrix.boundaryCoeffs()[sPatchI];

//        if (!ggiPatch_.master())
        {
            // This is master, interpolate shadow on this side
            const scalarField shadowInterpolatedIC = ggiPatch_.interpolate(shadowIC);
            const scalarField shadowInterpolatedBC = ggiPatch_.interpolate(shadowBC);

            // Set new internal coeffs using the data from the other side and
            // taking into account partially and fully uncovered faces
            const scalarField origIC = thisIC;
            thisIC = shadowInterpolatedBC;
            ggiPatch_.scaleForPartialCoverage(origIC, thisIC);
            ggiPatch_.scaleForPartialCoverage(origIC, thisIC);

            // Set new boundary coeffs using the data from the other side and
            // taking into account partially and fully uncovered faces
            const scalarField origBC = thisBC;
            thisBC = shadowInterpolatedIC;
            ggiPatch_.scaleForPartialCoverage(origBC, thisBC);
            ggiPatch_.scaleForPartialCoverage(origBC, thisBC);

            Info<< "This new internal coeffs: " << thisIC << nl
                << "This old internal coeffs: " << origIC << nl
                << "This shadow boundary coeffs" << shadowBC << nl << endl;
            
            Info<< "This new boundary coeffs: " << thisBC << nl
                << "This old boundary coeffs: " << origBC << nl
                << "This shadow internale coeffs" << shadowIC << nl << endl;

        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
