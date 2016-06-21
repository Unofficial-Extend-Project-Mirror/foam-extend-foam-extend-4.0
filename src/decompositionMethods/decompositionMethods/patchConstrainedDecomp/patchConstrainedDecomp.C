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

Description
    Decomposition given a cell-to-processor association in a file

\*---------------------------------------------------------------------------*/

#include "patchConstrainedDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchConstrainedDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        patchConstrainedDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchConstrainedDecomp::patchConstrainedDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh),
    dict_
    (
        decompositionDict.subDict
        (
            typeName + "Coeffs"
        )
    ),
    baseDecompPtr_
    (
        decompositionMethod::New(dict_, mesh)
    ),
    patchConstraints_(dict_.lookup("patchConstraints"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::patchConstrainedDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    labelList finalDecomp = baseDecompPtr_->decompose(points, pointWeights);

    // Impose the decomposition along patches
    forAll (patchConstraints_, i)
    {
        const label patchID =
            mesh_.boundaryMesh().findPatchID(patchConstraints_[i].first());

        const label procID = patchConstraints_[i].second();

        if (patchID < 0 || procID < 0 || procID > nProcessors_ - 1)
        {
            FatalErrorIn
            (
                "labelList patchConstrainedDecomp::decompose\n"
                "(\n"
                "    const pointField& points,\n"
                "    const scalarField& pointWeights\n"
                ")"
            )   << "Incorrect patch constraint definition for "
                << patchConstraints_[i]
                << abort(FatalError);
        }

        Info<< "Putting patch named " << patchConstraints_[i].first()
            << " index " << patchID
            << " to processor " << procID
            << endl;

        const labelList fc = mesh_.boundaryMesh()[patchID].faceCells();

        forAll (fc, fcI)
        {
            finalDecomp[fc[fcI]] = procID;
        }
    }

    fixCyclics(mesh_, finalDecomp);

    return finalDecomp;
}


// ************************************************************************* //
