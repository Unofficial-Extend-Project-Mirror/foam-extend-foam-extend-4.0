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

#include "cyclicFaPatch.H"
#include "coupledPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicFaPatch, 0);
addToRunTimeSelectionTable(faPatch, cyclicFaPatch, dictionary);

const scalar cyclicFaPatch::matchTol_
(
    debug::tolerances("patchFaceMatchTol", 1e-3)
);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicFaPatch::calcTransforms()
{
    if (size() > 0)
    {
        pointField half0Ctrs(size()/2);
        pointField half1Ctrs(size()/2);
        for (label i=0; i<size()/2; i++)
        {
            half0Ctrs[i] = this->edgeCentres()[i];
            half1Ctrs[i] = this->edgeCentres()[i+size()/2];
        }

        vectorField half0Normals(size()/2);
        vectorField half1Normals(size()/2);
 
        vectorField eN = edgeNormals()*magEdgeLengths();

        scalar maxMatchError = 0;
        label errorEdge = -1;

        for (label edgei = 0; edgei < size()/2; edgei++)
        {
            half0Normals[edgei] = eN[edgei];
            label nbrEdgei = edgei + size()/2;
            half1Normals[edgei] = eN[nbrEdgei];

            scalar magLe = mag(half0Normals[edgei]);
            scalar nbrMagLe = mag(half1Normals[edgei]);
            scalar avLe = (magLe + nbrMagLe)/2.0;

            if (magLe < ROOTVSMALL && nbrMagLe < ROOTVSMALL)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check. (note use of sqrt(VSMALL) since that is how mag
                // scales)
                half0Normals[edgei] = point(1, 0, 0);
                half1Normals[edgei] = half0Normals[edgei];
            }
            else if
            (
                mag(magLe - nbrMagLe)/avLe
              > matchTol_
            )
            {
                // Error in area matching.  Find largest error
                maxMatchError =
                    Foam::max(maxMatchError, mag(magLe - nbrMagLe)/avLe);
                errorEdge = edgei;
            }
            else
            {
                half0Normals[edgei] /= magLe;
                half1Normals[edgei] /= nbrMagLe;
            }
        }

        // Check for error in edge matching
        if (maxMatchError > matchTol_)
        {
            label nbrEdgei = errorEdge + size()/2;
            scalar magLe = mag(half0Normals[errorEdge]);
            scalar nbrMagLe = mag(half1Normals[errorEdge]);
            scalar avLe = (magLe + nbrMagLe)/2.0;

            FatalErrorIn
            (
                "cyclicFaPatch::calcTransforms()"
            )   << "edge " << errorEdge
                << " area does not match neighbour "
                << nbrEdgei << " by "
                << 100*mag(magLe - nbrMagLe)/avLe
                << "% -- possible edge ordering problem." << endl
                << "patch:" << name()
                << " my area:" << magLe
                << " neighbour area:" << nbrMagLe
                << " matching tolerance:" << matchTol_
                << endl
                << "Mesh edge:" << start() + errorEdge
                << endl
                << "Neighbour edge:" << start() + nbrEdgei
                << endl
                << "Other errors also exist, only the largest is reported. "
                << "Please rerun with cyclic debug flag set"
                << " for more information." << exit(FatalError);
        }

        // Calculate transformation tensors
        calcTransformTensors
        (
            half0Ctrs,
            half1Ctrs,
            half0Normals,
            half1Normals
        );

        // Check transformation tensors
        if (!parallel())
        {
            if (forwardT().size() > 1 || reverseT().size() > 1)
            {
                SeriousErrorIn
                (
                    "void cyclicFaPatch::calcTransforms()"
                )   << "Transformation tensor is not constant for the cyclic "
                    << "patch.  Please reconsider your setup and definition of "
                    << "cyclic boundaries." << endl;
            }
        }
    }
}


// Make patch weighting factors
void cyclicFaPatch::makeWeights(scalarField& w) const
{
    const scalarField& magL = magEdgeLengths();

    scalarField deltas = edgeNormals() & faPatch::delta();
    label sizeby2 = deltas.size()/2;

    scalar maxMatchError = 0;
    label errorEdge = -1;

    for (label edgei = 0; edgei < sizeby2; edgei++)
    {
        scalar avL = (magL[edgei] + magL[edgei + sizeby2])/2.0;

        if
        (
            mag(magL[edgei] - magL[edgei + sizeby2])/avL
          > matchTol_
        )
        {
            // Found error.  Look for largest matching error
            maxMatchError =
                Foam::max
                (
                    maxMatchError,
                    mag(magL[edgei] - magL[edgei + sizeby2])/avL
                );

            errorEdge = edgei;
        }

        scalar di = deltas[edgei];
        scalar dni = deltas[edgei + sizeby2];

        w[edgei] = dni/(di + dni);
        w[edgei + sizeby2] = 1 - w[edgei];
    }

    // Check for error in matching
    if (maxMatchError > polyPatch::matchTol_)
    {
        scalar avL = (magL[errorEdge] + magL[errorEdge + sizeby2])/2.0;

        FatalErrorIn("cyclicFaPatch::makeWeights(scalarField& w) const")
            << "edge " << errorEdge << " and " << errorEdge + sizeby2
            <<  " areas do not match by "
            << 100*mag(magL[errorEdge] - magL[errorEdge + sizeby2])/avL
            << "% -- possible edge ordering problem." << nl
            << "Cyclic area match tolerance = "
            << polyPatch::matchTol_ << " patch: " << name()
            << abort(FatalError);
    }
}


// Make patch edge - neighbour cell distances
void cyclicFaPatch::makeDeltaCoeffs(scalarField& dc) const
{
    scalarField deltas = edgeNormals() & faPatch::delta();
    label sizeby2 = deltas.size()/2;

    for (label edgei = 0; edgei < sizeby2; edgei++)
    {
        scalar di = deltas[edgei];
        scalar dni = deltas[edgei + sizeby2];

        dc[edgei] = 1.0/(di + dni);
        dc[edgei + sizeby2] = dc[edgei];
    }
}


void Foam::cyclicFaPatch::initGeometry()
{
    faPatch::initGeometry();
}


void Foam::cyclicFaPatch::calcGeometry()
{
    faPatch::calcGeometry();
    calcTransforms();
}


void Foam::cyclicFaPatch::initMovePoints(const pointField& p)
{
    faPatch::initMovePoints(p);
}


void Foam::cyclicFaPatch::movePoints(const pointField& p)
{
    faPatch::movePoints(p);
    calcTransforms();
}

// Return delta (P to N) vectors across coupled patch
tmp<vectorField> cyclicFaPatch::delta() const
{
    vectorField patchD = faPatch::delta();
    label sizeby2 = patchD.size()/2;

    tmp<vectorField> tpdv(new vectorField(patchD.size()));
    vectorField& pdv = tpdv();

    // Do the transformation if necessary
    if (parallel())
    {
        for (label edgei = 0; edgei < sizeby2; edgei++)
        {
            vector ddi = patchD[edgei];
            vector dni = patchD[edgei + sizeby2];

            pdv[edgei] = ddi - dni;
            pdv[edgei + sizeby2] = -pdv[edgei];
        }
    }
    else
    {
        for (label edgei = 0; edgei < sizeby2; edgei++)
        {
            vector ddi = patchD[edgei];
            vector dni = patchD[edgei + sizeby2];

            pdv[edgei] = ddi - transform(forwardT()[0], dni);
            pdv[edgei + sizeby2] = -transform(reverseT()[0], pdv[edgei]);
        }
    }

    return tpdv;
}


tmp<labelField> cyclicFaPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


tmp<labelField> cyclicFaPatch::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    tmp<labelField> tpnf(new labelField(this->size()));
    labelField& pnf = tpnf();

    label sizeby2 = this->size()/2;

    for (label edgei=0; edgei<sizeby2; edgei++)
    {
        pnf[edgei] = interfaceData[edgei + sizeby2];
        pnf[edgei + sizeby2] = interfaceData[edgei];
    }

    return tpnf;
}


tmp<labelField> cyclicFaPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    const unallocLabelList& edgeCells = this->faceCells();

    tmp<labelField> tpnf(new labelField(this->size()));
    labelField& pnf = tpnf();

    label sizeby2 = this->size()/2;

    for (label edgei=0; edgei<sizeby2; edgei++)
    {
        pnf[edgei] = iF[edgeCells[edgei + sizeby2]];
        pnf[edgei + sizeby2] = iF[edgeCells[edgei]];
    }

    return tpnf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
