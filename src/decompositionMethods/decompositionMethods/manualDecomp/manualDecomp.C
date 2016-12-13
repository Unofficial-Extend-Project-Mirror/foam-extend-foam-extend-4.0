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

#include "manualDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manualDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        manualDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manualDecomp::manualDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh),
    decompDataFile_
    (
        decompositionDict.subDict
        (
            word(decompositionDict.lookup("method"))
          + "Coeffs"
        ).lookup("dataFile")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::manualDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    labelIOList finalDecomp
    (
        IOobject
        (
            decompDataFile_,
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );

    // Check if the final decomposition is OK

    if (finalDecomp.size() != points.size())
    {
        FatalErrorIn
        (
            "manualDecomp::decompose(const pointField&, const scalarField&)"
        )   << "Size of decomposition list does not correspond "
            << "to the number of points.  Size: "
            << finalDecomp.size() << " Number of points: "
            << points.size()
            << ".\n" << "Manual decomposition data read from file "
            << decompDataFile_ << "." << endl
            << exit(FatalError);
    }

    if (min(finalDecomp) < 0 || max(finalDecomp) > nProcessors_ - 1)
    {
        FatalErrorIn
        (
            "manualDecomp::decompose(const pointField&, const scalarField&)"
        )   << "According to the decomposition, cells assigned to "
            << "impossible processor numbers.  Min processor = "
            << min(finalDecomp) << " Max processor = " << max(finalDecomp)
            << ".\n" << "Manual decomposition data read from file "
            << decompDataFile_ << "." << endl
            << exit(FatalError);
    }

    return finalDecomp;
}


// ************************************************************************* //
