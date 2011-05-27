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

Description
    Tetrahedral Finite Element matrix boundary treatment tools

\*---------------------------------------------------------------------------*/

#include "tetFemMatrices.H"
#include "tetPointRef.H"
#include "constraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Add boundary source for gradient-type conditions
template<class Type>
void tetFemMatrix<Type>::addBoundarySourceDiag()
{
    // Loop through all boundaries and fix the condition
    const FieldField<tetPolyPatchField, Type>& patches =
        psi().boundaryField();

    forAll (patches, patchI)
    {
        patches[patchI].addBoundarySourceDiag(*this);
    }
}


template<class Type>
void tetFemMatrix<Type>::storeBoundaryCoeffs() const
{
    // Make a map of all possible equations and only fix the ones that are
    // requested. Once the list of requested fixes is complete, collapse the
    // list. Note the use of list of pointers for the first part of operation
    // and creation of collapsed list by copying the pointers into a pointer
    // list.   HJ, 28/Feb/2001

    if (!boundaryConditionsSet_)
    {
        boundaryConditionsSet_ = true;

        // Loop through all boundaries and fix the condition
        const FieldField<tetPolyPatchField, Type>& patches =
            psi().boundaryField();

        forAll (patches, patchI)
        {
            patches[patchI].setBoundaryCondition(fixedEqns_);
        }

        // Loop through all fixed equations and grab the matrix
        labelList toc = fixedEqns_.toc();

        forAll (toc, eqnI)
        {
            fixedEqns_[toc[eqnI]].setMatrix(*this);
        }
    }
}


template<class Type>
void tetFemMatrix<Type>::setComponentBoundaryConditions
(
    const direction d,
    scalarField& psiCmpt,
    scalarField& sourceCmpt
)
{
    if (!boundaryConditionsSet_)
    {
        FatalErrorIn
        (
            "void tetFemMatrix<Type>::setComponentBoundaryConditions("
            "const direction& d, scalarField& psiCmpt, "
            "scalarField& sourceCmpt)"
        )   << "cannot reconstruct matrix: boundary conditions not set"
            << abort(FatalError);
    }

    // Loop through all fixed equations and set boundary condition
    labelList toc = fixedEqns_.toc();

    forAll (toc, eqnI)
    {
        fixedEqns_[toc[eqnI]].eliminateEquation(*this, d, sourceCmpt);
    }

    forAll (toc, eqnI)
    {
        fixedEqns_[toc[eqnI]].setSourceDiag(*this, d, psiCmpt, sourceCmpt);
    }
}


template<class Type>
void tetFemMatrix<Type>::reconstructMatrix()
{
    if (!boundaryConditionsSet_)
    {
        FatalErrorIn
        (
            "void tetFemMatrix<Type>::reconstructMatrix()"
        )   << "cannot reconstruct matrix: boundary conditions not set"
            << abort(FatalError);
    }

    // Loop through all fixed equations and reconstruct the matrix
    labelList toc = fixedEqns_.toc();

    forAll (toc, eqnI)
    {
        fixedEqns_[toc[eqnI]].reconstructMatrix(*this);
    }
}


template<class Type>
void tetFemMatrix<Type>::addCouplingCoeffs()
{
    // Diagonal

    if (hasDiag())
    {
        // Initialise diagonal transfer
        forAll (psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.initAddDiag(diag());
            }
        }

        // Add the diagonal
        forAll (psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.addDiag(diag());
            }
        }
    }

    // Upper triangle

    if (hasUpper())
    {
        // Initialise upper transfer
        forAll (psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.initAddUpperLower(upper());
            }
        }

        // Add the upper
        forAll (psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.addUpperLower(upper());
            }
        }
    }

    // Lower triangle

    if (hasLower())
    {
        // Initialise lower transfer
        forAll (psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.initAddUpperLower(lower());
            }
        }

        // Add the lower
        forAll (psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.addUpperLower(lower());
            }
        }
    }
}


template<class Type>
void tetFemMatrix<Type>::addCouplingSource(scalarField& sourceCmpt) const
{
    // Initialise source transfer
    forAll (psi_.boundaryField(), patchI)
    {
        const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.coupled())
        {
            ptf.initAddSource(sourceCmpt);
        }
    }

    // Add the source
    forAll (psi_.boundaryField(), patchI)
    {
        const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.coupled())
        {
            ptf.addSource(sourceCmpt);
        }
    }
}


template<class Type>
void tetFemMatrix<Type>::eliminateCouplingCoeffs()
{
    // Upper triangle
    if (hasUpper())
    {
        // Eliminate the upper
        forAll (psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.eliminateUpperLower(upper());
            }
        }
    }

    // Lower triangle

    if (hasLower())
    {
        // Eliminate the lower
        forAll (psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.eliminateUpperLower(lower());
            }
        }
    }
}


template<class Type>
tmp<Field<Type> > tetFemMatrix<Type>::distributeSource
(
    const Field<Type>& ef
) const
{
    const tetPolyMesh& mesh = psi().mesh();

    // Make a field over all points
    tmp<Field<Type> > tdistSu(mesh.nPoints(), pTraits<Type>::zero);
    Field<Type>& distSu = tdistSu();

    // Go through all the tets and add a quarter of the volume times ef
    // into each of the vertices. The cell integrals are prepared by the mesh
    labelList localToGlobalBuffer(mesh.maxNPointsForCell());
    labelList globalToLocalBuffer(lduAddr.size(), -1);

    scalarField coeffsBuffer
    (
        this->mesh.maxNTetsForCell()*tetPointRef::nVertices
    );

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        label nCellPoints =
            mesh.addressing
            (
                cellI,
                localToGlobalBuffer,
                globalToLocalBuffer
            );

        mesh.volIntegral(cellI, coeffsBuffer, globalToLocalBuffer);

        const Type& curEf = ef[cellI];

        for (label localI = 0; localI < nCellPoints; localI++)
        {
            label globalI = localToGlobalBuffer[localI];

            // Insert the source
            distSu[globalI] += coeffsBuffer[localI]*curEf;
        }

        // Clear addressing for element
        mesh.clearAddressing
        (
            cellI,
            nCellPoints,
            localToGlobalBuffer,
            globalToLocalBuffer
        );
    }

    return tdistSu;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
