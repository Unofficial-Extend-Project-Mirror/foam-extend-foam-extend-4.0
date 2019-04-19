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

#include "pointPatchInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyFvPatch.H"
#include "demandDrivenData.H"
#include "coupledPointPatchFields.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointPatchInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void pointPatchInterpolation::makePatchPatchAddressing() const
{
    if (debug)
    {
        Info<< "pointPatchInterpolation::makePatchPatchAddressing() : "
            << "constructing boundary addressing"
            << endl;
    }

    if
    (
        !patchInterpolators_.empty()
      || patchPatchPointsPtr_
      || patchPatchPointConstraintPointsPtr_
      || patchPatchPointConstraintTensorsPtr_
    )
    {
        FatalErrorIn
        (
            "void pointPatchInterpolation::makePatchPatchAddressing() const"
        )   << "Patch-patch point interpolation data already calculated."
            << abort(FatalError);
    }

    const fvBoundaryMesh& bm = fvMesh_.boundary();
    const pointBoundaryMesh& pbm = pointMesh::New(fvMesh_).boundary();

    // first count the total number of patch-patch points

    label nPatchPatchPoints = 0;

    forAll (bm, patchi)
    {
        if(!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
        {
            nPatchPatchPoints += bm[patchi].patch().boundaryPoints().size();
        }
    }


    // Go through all patches and mark up the external edge points
    Map<label> patchPatchPointSet(2*nPatchPatchPoints);

    // Allocate patch patch points
    patchPatchPointsPtr_ = new labelList(nPatchPatchPoints);
    labelList& patchPatchPoints = *patchPatchPointsPtr_;

    List<pointConstraint> patchPatchPointConstraints(nPatchPatchPoints);

    label pppi = 0;

    forAll (bm, patchi)
    {
        if(!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
        {
            const labelList& bp = bm[patchi].patch().boundaryPoints();
            const labelList& meshPoints = bm[patchi].patch().meshPoints();

            forAll (bp, pointi)
            {
                label ppp = meshPoints[bp[pointi]];

                Map<label>::iterator iter = patchPatchPointSet.find(ppp);

                if (iter == patchPatchPointSet.end())
                {
                    patchPatchPointSet.insert(ppp, pppi);
                    patchPatchPoints[pppi] = ppp;

                    pbm[patchi].applyConstraint
                    (
                        bp[pointi],
                        patchPatchPointConstraints[pppi]
                    );
                    pppi++;
                }
                else
                {
                    pbm[patchi].applyConstraint
                    (
                        bp[pointi],
                        patchPatchPointConstraints[iter()]
                    );
                }
            }
        }
    }

    nPatchPatchPoints = pppi;
    patchPatchPoints.setSize(nPatchPatchPoints);
    patchPatchPointConstraints.setSize(nPatchPatchPoints);

    // Allocate patch patch point constraint data
    patchPatchPointConstraintPointsPtr_ = new labelList(nPatchPatchPoints);
    labelList& patchPatchPointConstraintPoints =
        *patchPatchPointConstraintPointsPtr_;

    patchPatchPointConstraintTensorsPtr_ = new tensorField(nPatchPatchPoints);
    tensorField& patchPatchPointConstraintTensors =
        *patchPatchPointConstraintTensorsPtr_;


    label nConstraints = 0;

    forAll (patchPatchPointConstraints, i)
    {
        if (patchPatchPointConstraints[i].first() != 0)
        {
            patchPatchPointConstraintPoints[nConstraints] =
                patchPatchPoints[i];

            patchPatchPointConstraintTensors[nConstraints] =
                patchPatchPointConstraints[i].constraintTransformation();

            nConstraints++;
        }
    }

    patchPatchPointConstraintPoints.setSize(nConstraints);
    patchPatchPointConstraintTensors.setSize(nConstraints);

    // Set patch interpolators
    patchInterpolators_.setSize(bm.size());

    forAll (bm, patchi)
    {
        patchInterpolators_.set
        (
            patchi,
            new primitivePatchInterpolation(bm[patchi].patch())
        );
    }

    if (debug)
    {
        Info<< "pointPatchInterpolation::makePatchPatchAddressing() : "
            << "finished constructing boundary addressing"
            << endl;
    }
}


void pointPatchInterpolation::makePatchPatchWeights() const
{
    if (debug)
    {
        Info<< "pointPatchInterpolation::makePatchPatchWeights() : "
            << "constructing boundary weighting factors"
            << endl;
    }

    if (patchPatchPointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void pointPatchInterpolation::makePatchPatchWeights() const"
        )   << "Patch-patch point interpolation weights already calculated."
            << abort(FatalError);
    }

    // Get patch patch points
    const labelList& patchPatchPts = patchPatchPoints();

    // Allocate weights
    patchPatchPointWeightsPtr_ = new scalarListList(patchPatchPts.size());
    scalarListList& patchPatchPointWeights = *patchPatchPointWeightsPtr_;

    // Note: Do not use fvMesh functionality, because of update order
    // The polyMesh is already up-to-date - use this instead
    // HJ, 18/Feb/2011

    const labelListList& pf = fvMesh_.pointFaces();
    const fvBoundaryMesh& bm = fvMesh_.boundary();

    pointScalarField sumWeights
    (
        IOobject
        (
            "sumWeights",
            fvMesh_.polyMesh::instance(),
            fvMesh_
        ),
        pointMesh::New(fvMesh_),
        dimensionedScalar("zero", dimless, 0)
    );

    forAll (patchPatchPts, pointi)
    {
        const label curPoint = patchPatchPts[pointi];
        const labelList& curFaces = pf[curPoint];

        patchPatchPointWeights[pointi].setSize(curFaces.size());
        scalarList& pw = patchPatchPointWeights[pointi];

        label nFacesAroundPoint = 0;

        const vector& pointLoc = fvMesh_.points()[curPoint];

        forAll (curFaces, facei)
        {
            if (!fvMesh_.isInternalFace(curFaces[facei]))
            {
                label patchi =
                    fvMesh_.boundaryMesh().whichPatch(curFaces[facei]);

                if (!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
                {
                    pw[nFacesAroundPoint] =
                        1.0/mag
                        (
                            pointLoc
                            - fvMesh_.boundaryMesh()[patchi].faceCentres()
                                [bm[patchi].patch().whichFace(curFaces[facei])]
                        );

                    nFacesAroundPoint++;
                }
            }
        }

        // Reset the sizes of the local weights
        pw.setSize(nFacesAroundPoint);

        // Collect the sum of weights for parallel correction
        sumWeights[curPoint] += sum(pw);
    }

    // Do parallel correction of weights

    // Update coupled boundaries
    forAll (sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            sumWeights.boundaryField()[patchi].initAddField();
        }
    }

    forAll (sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            sumWeights.boundaryField()[patchi].addField
            (
                sumWeights.internalField()
            );
        }
    }

    // Re-scale the weights for the current point
    forAll (patchPatchPts, pointi)
    {
        scalarList& pw = patchPatchPointWeights[pointi];
        scalar sumw = sumWeights[patchPatchPts[pointi]];

        forAll (pw, facei)
        {
            pw[facei] /= sumw;
        }
    }

    if (debug)
    {
        Info<< "pointPatchInterpolation::makePatchPatchWeights() : "
            << "finished constructing boundary weighting factors"
            << endl;
    }
}


void pointPatchInterpolation::clearOut() const
{
    // Delete demand driven data
    deleteDemandDrivenData(patchPatchPointsPtr_);
    deleteDemandDrivenData(patchPatchPointConstraintPointsPtr_);
    deleteDemandDrivenData(patchPatchPointConstraintTensorsPtr_);
    deleteDemandDrivenData(patchPatchPointWeightsPtr_);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

pointPatchInterpolation::pointPatchInterpolation(const fvMesh& vm)
:
    fvMesh_(vm),
    patchInterpolators_(),
    patchPatchPointsPtr_(nullptr),
    patchPatchPointConstraintPointsPtr_(nullptr),
    patchPatchPointConstraintTensorsPtr_(nullptr),
    patchPatchPointWeightsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

pointPatchInterpolation::~pointPatchInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const PtrList<primitivePatchInterpolation>&
pointPatchInterpolation::patchInterpolators() const
{
    if (patchInterpolators_.empty())
    {
        makePatchPatchAddressing();
    }

    return patchInterpolators_;
}


const labelList& pointPatchInterpolation::patchPatchPoints() const
{
    if (!patchPatchPointsPtr_)
    {
        makePatchPatchAddressing();
    }

    return *patchPatchPointsPtr_;
}


const labelList&
pointPatchInterpolation::patchPatchPointConstraintPoints() const
{
    if (!patchPatchPointConstraintPointsPtr_)
    {
        makePatchPatchAddressing();
    }

    return *patchPatchPointConstraintPointsPtr_;
}


const tensorField&
pointPatchInterpolation::patchPatchPointConstraintTensors() const
{
    if (!patchPatchPointConstraintTensorsPtr_)
    {
        makePatchPatchAddressing();
    }

    return *patchPatchPointConstraintTensorsPtr_;
}


const scalarListList& pointPatchInterpolation::patchPatchPointWeights() const
{
    if (!patchPatchPointWeightsPtr_)
    {
        makePatchPatchWeights();
    }

    return *patchPatchPointWeightsPtr_;
}


bool pointPatchInterpolation::movePoints()
{
    // Clear out patch interpolators to force recalculation
    patchInterpolators_.clear();

    // Clear out demand driven data
    clearOut();

    return true;
}


void pointPatchInterpolation::updateMesh()
{
    // Clear out patch interpolators to force recalculation
    patchInterpolators_.clear();

    // Clear out demand driven data
    clearOut();
}


// Specialisation of applyCornerConstraints for scalars because
// no constraint need be applied
template<>
void pointPatchInterpolation::applyCornerConstraints<scalar>
(
    GeometricField<scalar, pointPatchField, pointMesh>& pf
) const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
