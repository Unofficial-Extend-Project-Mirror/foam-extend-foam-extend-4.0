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

\*---------------------------------------------------------------------------*/

#include "oversetFvPatchField.H"
#include "oversetFvPatch.H"
#include "fvMatrix.H"
#include "demandDrivenData.H"

#include "IndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
void oversetFvPatchField<Type>::setCoupledFringe
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    const bool coupled
)
{
    forAll(psi.boundaryField(), patchI)
    {
        fvPatchField<Type>& psip = psi.boundaryField()[patchI];

        if (isA<oversetFvPatchField<Type> >(psip))
        {
            oversetFvPatchField<Type>& opf =
                refCast<oversetFvPatchField<Type> >(psip);

            opf.setCoupledFringe(coupled);
        }
    }
}


template<class Type>
void oversetFvPatchField<Type>::oversetInterpolate
(
    GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    // Loop through boundary field and find overset patch
    forAll(psi.boundaryField(), patchI)
    {
        fvPatchField<Type>& psip = psi.boundaryField()[patchI];

        if (isA<oversetFvPatchField<Type> >(psip))
        {
            oversetFvPatchField<Type>& opf =
                refCast<oversetFvPatchField<Type> >(psip);

            // Store old settings: coupledFringe and conservativeCorrection
            const bool cfOrig = opf.coupled();
            const bool ccOrig = opf.conservativeCorrection();

            // Switch on coupled fringe and switch off conservative correction
            // to allow manual interpolation
            opf.setCoupledFringe(true);
            opf.setConservativeCorrection(false);

            // Perform overset interpolation
            opf.initEvaluate(Pstream::defaultComms()); // Performs interpolation
            opf.evaluate(Pstream::defaultComms()); // Sets hole cell values

            // Finally, revert to original settings
            opf.setCoupledFringe(cfOrig);
            opf.setConservativeCorrection(ccOrig);
        }
    }

    // After we have performed overset interpolation, we need to make sure that
    // the data correct data is transferred across processor boundaries
    psi.boundaryField().evaluateCoupled();
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void oversetFvPatchField<Type>::setHoleValues(Field<Type2>& f) const
{
    const labelList& dc = oversetPatch_.overset().holeCells();

    forAll (dc, dcI)
    {
        f[dc[dcI]] = holeCellValue_;
    }
}


template<class Type>
template<class Type2>
void oversetFvPatchField<Type>::setAcceptorValues(Field<Type2>& f) const
{
    // Get acceptor values by interpolation. Note that we assume that the field
    // we are operating on uses the same overset interpolation scheme as the
    // field referred by this fvPatchField. This might be problematic for GFM
    // interpolation schemes becayse f can be some arbitrary field depending on
    // the linear solver that is used. VV, 15/Feb/2017.
    Field<Type2> accValues = oversetPatch_.overset().interpolate
    (
        f,
        this->dimensionedInternalField().name()
    );

    // Get acceptor addressing
    const labelList& acceptors = oversetPatch_.overset().acceptorCells();

    // Check sizes
    if (accValues.size() != acceptors.size())
    {
        FatalErrorIn
        (
            "oversetFvPatchField<Type>::"
            "setAcceptorValues(Field<Type2>& f) const"
        )   << "Bad sizes: " << accValues.size()
            << " and " << acceptors.size()
            << abort(FatalError);
    }

    forAll (acceptors, accI)
    {
        f[acceptors[accI]] = accValues[accI];
    }
}


template<class Type>
void oversetFvPatchField<Type>::correctDiag
(
    fvMatrix<Type>& eqn
) const
{
    // Get access to diagonal
    scalarField& diag = eqn.diag();
    Field<Type>& source = eqn.source();

    const labelList& holeCells = oversetPatch_.overset().holeCells();
    const labelList& acceptorCells = oversetPatch_.overset().acceptorCells();

    label nLiveCells = diag.size() - holeCells.size() - acceptorCells.size();

    reduce(nLiveCells, sumOp<label>());

    // Estimate diagonal in live cells
    scalar liveDiag = 1;

    if (nLiveCells > 0)
    {
        liveDiag = gSumMag(diag)/nLiveCells;

        // Correct for sign
        liveDiag *= sign(gMax(diag));
    }
    else
    {
        FatalErrorIn
        (
            "void oversetFvPatchField<Type>::correctDiag\n"
            "(\n"
            "    fvMatrix<Type>& eqn\n"
            ") const"
        )   << "No live cells in matrix"
            << abort(FatalError);
    }

    // Fix diagonal if missing in hole cells
    forAll (holeCells, hcI)
    {
        if (mag(diag[holeCells[hcI]]) < SMALL)
        {
            diag[holeCells[hcI]] = liveDiag;
        }
    }

    // Fix diagonal (if missing) and source in acceptor cells
    forAll (acceptorCells, acI)
    {
        if (mag(diag[acceptorCells[acI]]) < SMALL)
        {
            diag[acceptorCells[acI]] = liveDiag;
        }

        source[acceptorCells[acI]] = pTraits<Type>::zero;
    }
}


template<class Type>
void oversetFvPatchField<Type>::correctOffDiag
(
    fvMatrix<Type>& eqn
) const
{
    // Kill off-diagonal coefficient in all hole and acceptor faces
    // Collect off-diagonal coefficients for all fringe faces

    // Get reference to overset mesh and fvMesh
    const oversetMesh& om = oversetPatch_.overset();
    const fvMesh& mesh = om.mesh();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    // 1 Hole internal faces
    const labelList& holeInternalFaces = om.holeInternalFaces();

    // 2 Acceptor internal faces
    const labelList& acceptorInternalFaces = om.acceptorInternalFaces();

    // 3 Hole faces
    const labelList& holeFaces = om.holeFaces();

    // 4 Fringe faces
    const labelList& fringeFaces = om.fringeFaces();
    const boolList& fringeFaceFlips = om.fringeFaceFlips();

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    if (eqn.symmetric())
    {
        scalarField& upper = eqn.upper();
        const label nInternalFaces = upper.size();

        // 1 Hole internal faces
        forAll (holeInternalFaces, hifI)
        {
            const label curFace = holeInternalFaces[hifI];

            if (curFace < nInternalFaces)
            {
                upper[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchI = boundaryMesh.whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchI].empty())
                {
                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(curFace);

                    eqn.internalCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                    eqn.boundaryCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 2 Acceptor internal faces
        forAll (acceptorInternalFaces, hifI)
        {
            const label curFace = acceptorInternalFaces[hifI];

            if (curFace < nInternalFaces)
            {
                upper[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                const label patchI = boundaryMesh.whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchI].empty())
                {
                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(curFace);

                    eqn.internalCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                    eqn.boundaryCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 3 Hole faces
        forAll (holeFaces, faceI)
        {
            const label curFace = holeFaces[faceI];

            if (curFace < nInternalFaces)
            {
                upper[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                const label patchI = boundaryMesh.whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchI].empty())
                {
                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(curFace);

                    eqn.internalCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                    eqn.boundaryCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 4 Fringe faces
        fringeUpperCoeffs_.setSize(fringeFaces.size(), 0);
        fringeLowerCoeffs_.setSize(fringeFaces.size(), 0);

        forAll (fringeFaces, fringeI)
        {
            const label curFace = fringeFaces[fringeI];

            if (curFace < nInternalFaces)
            {
                // Internal face: kill the off-diagonal coefficient in
                // the live cell

                // Since the matrix is symmetric, there is no need to
                // distinguish between lower and upper coefficient
                fringeUpperCoeffs_[fringeI] = upper[curFace];
                fringeLowerCoeffs_[fringeI] = upper[curFace];

                // Kill original coefficient only if coupledFringe is used. If
                // the coupledFringe_ = false (explicit fringe update), we need
                // to have the matrix coefficients in order to use
                // fvMatrix::setValues member function.  VV, 23/Sep/2016.
                if (coupledFringe_)
                {
                    upper[curFace] = 0;
                }
            }
            else
            {
                // Coupled boundary
                const label patchI = boundaryMesh.whichPatch(curFace);

                // First check if the patch is coupled
                if (psi.boundaryField()[patchI].coupled())
                {
                    // Copy boundary/internal coeffs into fringe upper/lower

                    // Note: on coupled boundaries, all coefficients are
                    // identical. We can take the first component
                    // HJ, 30/May/2013

                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(curFace);

                    fringeUpperCoeffs_[fringeI] =
                        component
                        (
                            eqn.boundaryCoeffs()[patchI][patchFaceI],
                            0
                        );
                    fringeLowerCoeffs_[fringeI] =
                        component
                        (
                            eqn.internalCoeffs()[patchI][patchFaceI],
                            0
                        );

                    // Check if the acceptor is on this side
                    // 1. If it is: kill corresponding internal/boundary
                    //    coeffs
                    // 2. If it isn't: do nothing, this is a live cell
                    if (fringeFaceFlips[fringeI])
                    {
                        // Kill coeffs only if coupled fringe is used
                        if (coupledFringe_)
                        {
                            eqn.internalCoeffs()[patchI][patchFaceI] =
                                pTraits<Type>::zero;
                            eqn.boundaryCoeffs()[patchI][patchFaceI] =
                                pTraits<Type>::zero;
                        }
                    }
                }
            }
        }
    }
    else if (eqn.asymmetric())
    {
        scalarField& upper = eqn.upper();
        scalarField& lower = eqn.lower();
        const label nInternalFaces = upper.size();

        // 1 Hole internal faces
        forAll (holeInternalFaces, hifI)
        {
            const label curFace = holeInternalFaces[hifI];

            if (curFace < nInternalFaces)
            {
                upper[curFace] = 0;
                lower[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                const label patchI = boundaryMesh.whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchI].empty())
                {
                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(curFace);

                    eqn.internalCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                    eqn.boundaryCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 2 Acceptor internal faces
        forAll (acceptorInternalFaces, hifI)
        {
            const label curFace = acceptorInternalFaces[hifI];

            if (curFace < nInternalFaces)
            {
                upper[curFace] = 0;
                lower[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                const label patchI = boundaryMesh.whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchI].empty())
                {
                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(curFace);

                    eqn.internalCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                    eqn.boundaryCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 3 Hole faces
        forAll (holeFaces, faceI)
        {
            const label curFace = holeFaces[faceI];

            if (curFace < nInternalFaces)
            {
                upper[curFace] = 0;
                lower[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                const label patchI = boundaryMesh.whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchI].empty())
                {
                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(curFace);

                    eqn.internalCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                    eqn.boundaryCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 4 Fringe faces
        fringeUpperCoeffs_.setSize(fringeFaces.size(), 0);
        fringeLowerCoeffs_.setSize(fringeFaces.size(), 0);

        forAll (fringeFaces, fringeI)
        {
            const label curFace = fringeFaces[fringeI];

            if (curFace < nInternalFaces)
            {
                // Internal face: kill the off-diagonal coefficient in
                // the live cell

                // Since the matrix is asymmetric, lower and upper
                // coefficients need to be distinguished
                fringeUpperCoeffs_[fringeI] = upper[curFace];
                fringeLowerCoeffs_[fringeI] = lower[curFace];

                // Kill original coefficient only if coupledFringe is used. If
                // the coupledFringe_ = false (explicit fringe update), we need
                // to have the matrix coefficients in order to use
                // fvMatrix::setValues member function.  VV, 23/Sep/2016.
                if (coupledFringe_)
                {
                    upper[curFace] = 0;
                    lower[curFace] = 0;
                }
            }
            else
            {
                // Coupled boundary
                const label patchI = boundaryMesh.whichPatch(curFace);

                // First check if the patch is coupled
                if (psi.boundaryField()[patchI].coupled())
                {
                    // Copy boundary/internal coeffs into fringe upper/lower

                    // Note: on coupled boundaries, all coefficients are
                    // identical. We can take the first component
                    // HJ, 30/May/2013

                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(curFace);

                    fringeUpperCoeffs_[fringeI] =
                        component
                        (
                            eqn.boundaryCoeffs()[patchI][patchFaceI],
                            0
                        );
                    fringeLowerCoeffs_[fringeI] =
                        component
                        (
                            eqn.internalCoeffs()[patchI][patchFaceI],
                            0
                        );

                    // Check if the acceptor is on this side
                    // 1. If it is: kill corresponding internal/boundary
                    //    coeffs
                    // 2. If it isn't: do nothing, this is a live cell
                    if (fringeFaceFlips[fringeI])
                    {
                        // Kill coeffs only if coupled fringe is used
                        if (coupledFringe_)
                        {
                            eqn.internalCoeffs()[patchI][patchFaceI] =
                                pTraits<Type>::zero;
                            eqn.boundaryCoeffs()[patchI][patchFaceI] =
                                pTraits<Type>::zero;
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
void oversetFvPatchField<Type>::correctFringeConservation
(
    scalarField& psiInterpolated,
    const lduMatrix& m
) const
{
    // Get necessary references
    const oversetMesh& om = oversetPatch_.overset();
    const polyBoundaryMesh& boundaryMesh = om.mesh().boundaryMesh();

    const labelList& fringeFaces = om.fringeFaces();
    const boolList& fringeFaceFlips = om.fringeFaceFlips();

    // Owner/neighbour and fringe addressing
    const unallocLabelList& own = m.lduAddr().lowerAddr();
    const unallocLabelList& nei = m.lduAddr().upperAddr();
    const label nInternalFaces = nei.size();

    // Initialise variables for fringe continuity error and sum of off diagonal
    // elements for fringe faces
    scalar fringeConservationError = 0;
    scalar fringeSumOffDiag = 0;

    // Loop through fringe faces and calculate flux balance
    forAll (fringeFaces, fringeI)
    {
        // Get current face index
        const label& curFace = fringeFaces[fringeI];
        const bool& curFlip = fringeFaceFlips[fringeI];

        if (curFace < nInternalFaces)
        {
            // Get fringe upper/lower coeffs
            const scalar& ufc = fringeUpperCoeffs_[fringeI];
            const scalar& lfc = fringeLowerCoeffs_[fringeI];

            // Calculate current flux
            const scalar curFlux =
                ufc*psiInterpolated[nei[curFace]]
              - lfc*psiInterpolated[own[curFace]];

            if (curFlip)
            {
                // Face points into a live cell: acceptor is owner. Flip the
                // sign of the flux and add lower coefficient to sum off diag
                fringeConservationError -= curFlux;
                fringeSumOffDiag += lfc;
            }
            else
            {
                // Face points into an acceptor: acceptor is neighbour. Flux is
                // correctly oriented, add upper coefficient to sum off diag
                fringeConservationError += curFlux;
                fringeSumOffDiag += ufc;
            }
        }
        else
        {
            // Coupled boundary. Note: each side (processor) takes care of its
            // own contribution which will later be combined
            const label patchI = boundaryMesh.whichPatch(curFace);
            const polyPatch& curPolyPatch = boundaryMesh[patchI];

            if (isA<processorPolyPatch>(curPolyPatch))
            {
                // Get processor polyPatch
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(curPolyPatch);

                // Get face cell index
                const label patchFaceI = curPolyPatch.whichFace(curFace);
                const label faceCellI = curPolyPatch.faceCells()[patchFaceI];

                // Get fringe lower coeffs
                const scalar& lfc = fringeLowerCoeffs_[fringeI];

                // Get interpolated psi value for this face cell
                const scalar& psiOwn = psiInterpolated[faceCellI];

                // Note: not entirely sure this is ok - what if lower/upper
                // coeffs for the same processor face on opposing sides are
                // different? Need to think about it a bit harder. Should work
                // for symmetric (pressure) matrices though. VV, 22/Jan/2016.
                // Note: reverse sign for faces on coupled boundaries
                if (curFlip)
                {
                    // Acceptor is on this side, add part of the flux from the
                    // acceptor on this side.
                    fringeConservationError -= lfc*psiOwn;

                    // Only owner processor patch contributes to sum off diag
                    if (procPatch.owner())
                    {
                        fringeSumOffDiag -= lfc;
                    }
                }
                else
                {
                    // Acceptor is on the other side, add part of the flux from
                    // the live cell on this side.
                    fringeConservationError += lfc*psiOwn;

                    // Only owner processor patch contributes to sum off diag
                    if (procPatch.owner())
                    {
                        fringeSumOffDiag -= lfc;
                    }
                }
            }
        }
    }

    // Perform global reduce on fringe conservation error and sum off diag
    reduce(fringeConservationError, sumOp<scalar>());
    reduce(fringeSumOffDiag, sumOp<scalar>());

    // Calculate bulk correction
    const scalar delta = -fringeConservationError/fringeSumOffDiag;

    // Get acceptors
    const labelList& acceptors = om.acceptorCells();

    // Correct interpolated acceptor values
    forAll (acceptors, accI)
    {
        psiInterpolated[acceptors[accI]] += delta;
    }
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    coupledFringe_(false),
    conservativeCorrection_(false),
    setHoleCellValue_(false),
    holeCellValue_(pTraits<Type>::zero),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    coupledFringe_(false),
    conservativeCorrection_(false),
    setHoleCellValue_(false),
    holeCellValue_(pTraits<Type>::zero),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    coupledFringe_(dict.lookup("coupledFringe")),
    conservativeCorrection_
    (
        dict.lookupOrDefault<Switch>("conservativeCorrection", false)
    ),
    setHoleCellValue_(dict.lookup("setHoleCellValue")),
    holeCellValue_(pTraits<Type>(dict.lookup("holeCellValue"))),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{
    if (!isType<oversetFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "oversetFvPatchField<Type>::oversetFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    coupledFringe_(ptf.coupledFringe_),
    conservativeCorrection_(ptf.conservativeCorrection_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{
    if (!isType<oversetFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "oversetFvPatchField<Type>::oversetFvPatchField\n"
            "(\n"
            "    const oversetFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf
)
:
    oversetLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    oversetPatch_(refCast<const oversetFvPatch>(ptf.patch())),
    coupledFringe_(ptf.coupledFringe_),
    conservativeCorrection_(ptf.conservativeCorrection_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    oversetPatch_(refCast<const oversetFvPatch>(ptf.patch())),
    coupledFringe_(ptf.coupledFringe_),
    conservativeCorrection_(ptf.conservativeCorrection_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{}


// * * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * //

template<class Type>
oversetFvPatchField<Type>::~oversetFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > oversetFvPatchField<Type>::patchNeighbourField() const
{
    // Warning: returning own patch field, as "neighbour" does not exist
    return *this;
}


template<class Type>
void oversetFvPatchField<Type>::updateCoeffs()
{
    // Fix the value in hole cells.  Probably unnecessary.  HJ, 25/Jun/2013
    if (setHoleCellValue_)
    {
        // Get non-const access to internal field
        Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

        this->setHoleValues(psiI);
    }

    // Clear fringe coefficients
    fringeUpperCoeffs_.clear();
    fringeLowerCoeffs_.clear();

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
void oversetFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Note that the acceptor cells should not be updated when:
    // 1. The fringe is not coupled (i.e. explicit update) because this update
    //    is called at the end of fvMatrix::solve and would mess up the equation
    //    fluxes (e.g. pEqn.flux() would not yield fully conservative flux)
    // 2. The fringe is coupled and the conservative correction had been used
    //    during the solution process. This would also mess up the equation
    //    fluxes.
    if (coupledFringe_ && !conservativeCorrection_)
    {
        if (oversetMesh::debug)
        {
            Info<< "Explicitly correcting fringe values for field: "
                << this->dimensionedInternalField().name() << endl;
        }

        // Get non-constant access to internal field
        Field<Type>& psi = const_cast<Field<Type>&>(this->internalField());

        // Set acceptor values
        this->setAcceptorValues(psi);
    }
}


template<class Type>
void oversetFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (setHoleCellValue_)
    {
        // Get non-constant access to internal field
        Field<Type>& psi = const_cast<Field<Type>&>(this->internalField());

        this->setHoleValues(psi);
    }
}


template<class Type>
void oversetFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& eqn
)
{
    // Eliminate unnecesary off-diagonal coefficients
    this->correctOffDiag(eqn);

    // Build matrix diagonal for cells where it is missing
    this->correctDiag(eqn);

    // Set values in hole cells
    const labelList& holeCells = oversetPatch_.overset().holeCells();

    // Correct equation for hole cells
    Field<Type> holeCellsPsi
    (
        holeCells.size(),
        holeCellValue_
    );

    eqn.setValues(holeCells, holeCellsPsi);

    // If the fringe is not coupled, set values in acceptor cells and do not
    // perform overset interpolation during the solution process. Note that the
    // fringeUpper/Lower coefficients are collected in correctOffDiag member
    // function, which is needed for correct flux reconstruction in explicit
    // fringe update (coupledFringe_ = false)
    if (!coupledFringe_)
    {
        Info<< "Setting values in acceptor cells before the solution for "
            << "field: " << eqn.psi().name() << endl;

        // Get acceptor values by interpolation
        Field<Type> accValues =
            oversetPatch_.overset().interpolate(eqn.psi(), eqn.psi().name());

        // Get acceptor addressing
        const labelList& accCells = oversetPatch_.overset().acceptorCells();

        eqn.setValues(accCells, accValues);
    }
}


template<class Type>
void oversetFvPatchField<Type>::transformCoupleField
(
    scalarField& f,
    const direction cmpt
) const
{}


template<class Type>
void oversetFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011
    if (coupledFringe_)
    {
        // Do multiplication on fringe faces

        // Get non-const access to psi
        scalarField& psi = const_cast<scalarField&>(psiInternal);

        // Set acceptor values
        this->setAcceptorValues(psi);

        // Correct for global fringe mass conservation
        if (conservativeCorrection_)
        {
            correctFringeConservation(psi, m);
        }

        const unallocLabelList& own = m.lduAddr().lowerAddr();
        const unallocLabelList& nei = m.lduAddr().upperAddr();

        const labelList& fringeFaces = oversetPatch_.overset().fringeFaces();
        const boolList& fringeFaceFlips =
            oversetPatch_.overset().fringeFaceFlips();

        if (switchToLhs)
        {
            forAll (fringeFaces, fringeI)
            {
                const label& curFace = fringeFaces[fringeI];

                // Multiplication is only done for internal fringe faces. Live
                // cells having fringe faces on coupled boundaries are correctly
                // taken into account via internal/boundary coeffs.
                // HJ, 1/May/2015
                if (curFace < own.size())
                {
                    // Get addressing
                    const label& o = own[curFace];
                    const label& n = nei[curFace];

                    // Note change of sign in multiplication, because
                    // fringe coefficients belong to A
                    // HJ, 22/May/2013
                    if (fringeFaceFlips[fringeI])
                    {
                        // Face pointing into live cell
                        // Add overset off-diagonal contribution to live cell
                        result[n] -= fringeLowerCoeffs_[fringeI]*psi[o];
                    }
                    else
                    {
                        // Face pointing out of live cell
                        // Add overset off-diagonal contribution to live cell
                        result[o] -= fringeUpperCoeffs_[fringeI]*psi[n];
                    }
                }
            }
        }
        else
        {
            forAll (fringeFaces, fringeI)
            {
                const label& curFace = fringeFaces[fringeI];

                // Multiplication is only done for internal fringe faces. Live
                // cells having fringe faces on coupled boundaries are correctly
                // taken into account via internal/boundary coeffs.
                // HJ, 1/May/2015
                if (curFace < own.size())
                {
                    // Get addressing
                    const label& o = own[curFace];
                    const label& n = nei[curFace];

                    // Note change of sign in multiplication, because
                    // fringe coefficients belong to A
                    // HJ, 22/May/2013
                    if (fringeFaceFlips[fringeI])
                    {
                        // Face pointing into of live cell
                        // Add overset off-diagonal contribution to live cell
                        result[n] += fringeLowerCoeffs_[fringeI]*psi[o];
                    }
                    else
                    {
                        // Face pointing out live cell
                        // Add overset off-diagonal contribution to live cell
                        result[o] += fringeUpperCoeffs_[fringeI]*psi[n];
                    }
                }
            }
        }

        // Do acceptor cells
        // Get diagonal
        const scalarField& diag = m.diag();

        const labelList& acceptorCells =
            oversetPatch_.overset().acceptorCells();

        forAll (acceptorCells, acI)
        {
            const label& curCell = acceptorCells[acI];

            result[curCell] = -diag[curCell]*psi[curCell];
        }
    }
    // else
    // Fringe is explicit (coupledFringe_ = false) and the acceptor cells are
    // treated as fixedValue which have been set via overset interpolation in
    // manipulateMatrix member function
}


template<class Type>
void oversetFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{}


template<class Type>
void oversetFvPatchField<Type>::patchFlux
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& pFlux,
    const fvMatrix<Type>& matrix
) const
{
    // Get necessary references
    const oversetMesh& om = oversetPatch_.overset();
    const polyBoundaryMesh& boundaryMesh = om.mesh().boundaryMesh();

    // Cast away constness of the fvMatrix (we need to change fringe face
    // coefficients for coupled boundaries)
    fvMatrix<Type>& eqn = const_cast<fvMatrix<Type>&>(matrix);

    // Get addressing
    const unallocLabelList& own = matrix.lduAddr().lowerAddr();
    const unallocLabelList& nei = matrix.lduAddr().upperAddr();
    const label nInternalFaces = nei.size();

    // Get fringe addressing
    const labelList& fringeFaces = om.fringeFaces();
    const boolList& fringeFaceFlips = om.fringeFaceFlips();

    // Get fields
    const GeometricField<Type, fvPatchField, volMesh>& psi = matrix.psi();
    const Field<Type>& psiIn = psi.internalField();
    Field<Type>& fluxIn = pFlux.internalField();

    // Note that fringe corrects coefficients on internal faces
    forAll (fringeFaces, fringeI)
    {
        // Get addressing
        const label& faceI = fringeFaces[fringeI];

        // Multiplication is only done for internal fringe faces.
        // Live/acceptors cells having fringe faces on coupled boundaries
        // are correctly taken into account via internal/boundary coeffs on
        // corresponding processor boundary(see below).
        if (faceI < nInternalFaces)
        {
            fluxIn[faceI] = fringeUpperCoeffs_[fringeI]*psiIn[nei[faceI]]
                - fringeLowerCoeffs_[fringeI]*psiIn[own[faceI]];
        }
        else
        {
            // Coupled boundary treatment:
            // Rebuild internal/boundary coefficients from fringe
            // lower/upper coefficients for consistent calculation of fringe
            // face fluxes on coupled boundaries. This implies that the
            // overset patch must come before processor patches in boundary
            // field - this check is performed in oversetMesh::oversetMesh.
            // VV, 25/Jan/2015.

            const label patchI = boundaryMesh.whichPatch(faceI);

            // First check if the patch is coupled
            if (psi.boundaryField()[patchI].coupled())
            {
                // Restore boundary/internal coeffs from fringe upper/lower
                // if the acceptor is on this side. If it isn't,
                // boundary/internal coefficients are correct - there's
                // nothing to do.
                if (fringeFaceFlips[fringeI])
                {
                    const label patchFaceI =
                        boundaryMesh[patchI].whichFace(faceI);

                    eqn.internalCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::one*fringeLowerCoeffs_[fringeI];
                    eqn.boundaryCoeffs()[patchI][patchFaceI] =
                        pTraits<Type>::one*fringeUpperCoeffs_[fringeI];
                }
            }
        }
    }
}


template<class Type>
void oversetFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeKeyword("coupledFringe") << coupledFringe_
        << token::END_STATEMENT << nl;
    os.writeKeyword("conservativeCorrection") << conservativeCorrection_
        << token::END_STATEMENT << nl;

    os.writeKeyword("setHoleCellValue")
        << setHoleCellValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("holeCellValue")
        << holeCellValue_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
