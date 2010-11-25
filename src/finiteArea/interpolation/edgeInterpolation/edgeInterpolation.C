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
    Face to edge interpolation scheme. Included in faMesh.

\*---------------------------------------------------------------------------*/

#include "faMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "demandDrivenData.H"
#include "faPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(edgeInterpolation, 0);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void edgeInterpolation::clearOut()
{
    deleteDemandDrivenData(lPN_);
    deleteDemandDrivenData(weightingFactors_);
    deleteDemandDrivenData(differenceFactors_);
    deleteDemandDrivenData(correctionVectors_);
    deleteDemandDrivenData(skewCorrectionVectors_);
//     deleteDemandDrivenData(leastSquarePvectors_);
//     deleteDemandDrivenData(leastSquareNvectors_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

edgeInterpolation::edgeInterpolation(const faMesh& fam)
:
    faSchemes(fam()),
    faSolution(fam()),
    faMesh_(fam),
    lPN_(NULL),
    weightingFactors_(NULL),
    differenceFactors_(NULL),
    orthogonal_(false),
    correctionVectors_(NULL),
    skew_(true),
    skewCorrectionVectors_(NULL)
//     leastSquarePvectors_(NULL),
//     leastSquareNvectors_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

edgeInterpolation::~edgeInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const edgeScalarField& edgeInterpolation::lPN() const
{
    if (!lPN_)
    {
        makeLPN();
    }

    return (*lPN_);
}


const edgeScalarField& edgeInterpolation::weights() const
{
    if (!weightingFactors_)
    {
        makeWeights();
    }

    return (*weightingFactors_);
}


const edgeScalarField& edgeInterpolation::deltaCoeffs() const
{
    if (!differenceFactors_)
    {
        makeDeltaCoeffs();
    }

    return (*differenceFactors_);
}


bool edgeInterpolation::orthogonal() const
{
    if (orthogonal_ == false && !correctionVectors_)
    {
        makeCorrectionVectors();
    }

    return orthogonal_;
}


const edgeVectorField& edgeInterpolation::correctionVectors() const
{
    if (orthogonal())
    {
        FatalErrorIn("edgeInterpolation::correctionVectors()")
            << "cannot return correctionVectors; mesh is orthogonal"
            << abort(FatalError);
    }

    return (*correctionVectors_);
}


bool edgeInterpolation::skew() const
{
    if (skew_ == true && !skewCorrectionVectors_)
    {
        makeSkewCorrectionVectors();
    }

    return skew_;
}


const edgeVectorField& edgeInterpolation::skewCorrectionVectors() const
{
    if (!skew())
    {
        FatalErrorIn("edgeInterpolation::skewCorrectionVectors()")
            << "cannot return skewCorrectionVectors; mesh is now skewed"
            << abort(FatalError);
    }

    return (*skewCorrectionVectors_);
}


// const edgeVectorField& edgeInterpolation::leastSquarePvectors() const
// {
//     if (!leastSquarePvectors_)
//     {
//         makeLeastSquareVectors();
//     }

//     return (*leastSquarePvectors_);
// }


// const edgeVectorField& edgeInterpolation::leastSquareNvectors() const
// {
//     if (!leastSquareNvectors_)
//     {
//         makeLeastSquareVectors();
//     }

//     return (*leastSquareNvectors_);
// }


// Do what is neccessary if the mesh has moved
bool edgeInterpolation::movePoints()
{
    deleteDemandDrivenData(lPN_);
    deleteDemandDrivenData(weightingFactors_);
    deleteDemandDrivenData(differenceFactors_);

    orthogonal_ = false;
    deleteDemandDrivenData(correctionVectors_);

    skew_ = true;
    deleteDemandDrivenData(skewCorrectionVectors_);

//     deleteDemandDrivenData(leastSquarePvectors_);
//     deleteDemandDrivenData(leastSquareNvectors_);

    return true;
}


void edgeInterpolation::makeLPN() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeLPN() : "
            << "Constructing geodesic distance between points P and N"
            << endl;
    }


    lPN_ = new edgeScalarField
    (
        IOobject
        (
            "lPN",
            faSolution::time().constant(),
            faSolution::db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimLength
    );
    edgeScalarField& lPN = *lPN_;


    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll(owner, edgeI)
    {
        vector curSkewCorrVec = vector::zero;

        if (skew())
        {
            curSkewCorrVec = skewCorrectionVectors()[edgeI];
        }

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - curSkewCorrVec
              - faceCentres[owner[edgeI]]
            );

        scalar lEN = 
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + curSkewCorrVec
            );

        lPN.internalField()[edgeI] = (lPE + lEN);
    }


    forAll(lPN.boundaryField(), patchI)
    {
        mesh().boundary()[patchI].makeDeltaCoeffs
        (
            lPN.boundaryField()[patchI]
        );

        lPN.boundaryField()[patchI] = 1.0/lPN.boundaryField()[patchI];
    }


    if (debug)
    {
        Info<< "edgeInterpolation::makeLPN() : "
            << "Finished constructing geodesic distance PN"
            << endl;
    }
}


void edgeInterpolation::makeWeights() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeWeights() : "
            << "Constructing weighting factors for edge interpolation"
            << endl;
    }


    weightingFactors_ = new edgeScalarField
    (
        IOobject
        (
            "weightingFactors",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimless
    );
    edgeScalarField& weightingFactors = *weightingFactors_;


    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll(owner, edgeI)
    {
        vector curSkewCorrVec = vector::zero;

        if (skew())
        {
            curSkewCorrVec = skewCorrectionVectors()[edgeI];
        }

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - curSkewCorrVec
              - faceCentres[owner[edgeI]]
            );

        scalar lEN = 
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + curSkewCorrVec
            );

        weightingFactors.internalField()[edgeI] = 
            lEN
            /(
                lPE
#               ifdef BAD_MESH_STABILISATION
              + VSMALL
#               endif
              + lEN
            );
    }

    forAll(mesh().boundary(), patchI)
    {
        mesh().boundary()[patchI].makeWeights
        (
            weightingFactors.boundaryField()[patchI]
        );
    }

    if (debug)
    {
        Info<< "edgeInterpolation::makeWeights() : "
            << "Finished constructing weighting factors for face interpolation"
            << endl;
    }
}


void edgeInterpolation::makeDeltaCoeffs() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeDeltaCoeffs() : "
            << "Constructing differencing factors array for edge gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    differenceFactors_ = new edgeScalarField
    (
        IOobject
        (
            "differenceFactors_",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimless/dimLength
    );
    edgeScalarField& DeltaCoeffs = *differenceFactors_;
    scalarField& dc = DeltaCoeffs.internalField();


    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();
    const edgeVectorField& lengths = mesh().Le();

    const edgeList& edges = mesh().edges();
    const pointField& points = mesh().points();


    forAll(owner, edgeI)
    {
        // Edge normal - area normal
        vector edgeNormal = lengths[edgeI]^edges[edgeI].vec(points);
        edgeNormal /= mag(edgeNormal);


        // Unit delta vector
        vector unitDelta =
            faceCentres[neighbour[edgeI]]
          - faceCentres[owner[edgeI]];

        unitDelta -=
            edgeNormal*(edgeNormal&unitDelta);

        unitDelta /= mag(unitDelta);


        // Calc PN arc length
        vector curSkewCorrVec = vector::zero;
        
        if (skew())
        {
            curSkewCorrVec = skewCorrectionVectors()[edgeI];
        }

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - curSkewCorrVec
              - faceCentres[owner[edgeI]]
            );

        scalar lEN = 
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + curSkewCorrVec
            );

        scalar lPN = lPE + lEN;


        // Edge normal - area tangent
        edgeNormal = lengths[edgeI]/mag(lengths[edgeI]);

        // Delta coefficient as per definition
//         dc[edgeI] = 1.0/(lPN*(unitDelta & edgeNormal));

        // Stabilised form for bad meshes.  HJ, 23/Jul/2009
        dc[edgeI] = 1.0/max((lPN*(unitDelta & edgeNormal)), 0.05*lPN);
    }


    forAll(DeltaCoeffs.boundaryField(), patchI)
    {
        mesh().boundary()[patchI].makeDeltaCoeffs
        (
            DeltaCoeffs.boundaryField()[patchI]
        );
    }
}


void edgeInterpolation::makeCorrectionVectors() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeCorrectionVectors() : "
            << "Constructing non-orthogonal correction vectors"
            << endl;
    }

    correctionVectors_ = new edgeVectorField
    (
        IOobject
        (
            "correctionVectors",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimless
    );
    edgeVectorField& CorrVecs = *correctionVectors_;

    // Set local references to mesh data
    const areaVectorField& faceCentres = mesh().areaCentres();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const edgeVectorField& lengths = mesh().Le();
    const edgeScalarField& magLengths = mesh().magLe();

    const edgeList& edges = mesh().edges();
    const pointField& points = mesh().points();

    scalarField deltaCoeffs(owner.size());

    forAll(owner, edgeI)
    {
        // Edge normal - area normal
        vector edgeNormal = lengths[edgeI] ^ edges[edgeI].vec(points);
        edgeNormal /= mag(edgeNormal);

        // Unit delta vector
        vector unitDelta =
            faceCentres[neighbour[edgeI]]
          - faceCentres[owner[edgeI]];

        unitDelta -= edgeNormal*(edgeNormal & unitDelta);

        unitDelta /= mag(unitDelta);

        // Edge normal - area tangent
        edgeNormal = lengths[edgeI]/magLengths[edgeI];

        // Delta coeffs
        deltaCoeffs[edgeI] = 1.0/(unitDelta & edgeNormal);

        // Edge correction vector
        CorrVecs.internalField()[edgeI] =
            edgeNormal
          - deltaCoeffs[edgeI]*unitDelta;
    }

    forAll(CorrVecs.boundaryField(), patchI)
    {
        faePatchVectorField& patchCorrVecs = CorrVecs.boundaryField()[patchI];

        patchCorrVecs = vector::zero;
    }

    scalar NonOrthogCoeff = 0.0;

    if (owner.size() > 0)
    {
        scalarField sinAlpha = deltaCoeffs*mag(CorrVecs.internalField());

        forAll(sinAlpha, edgeI)
        {
            sinAlpha[edgeI] = max(-1, min(sinAlpha[edgeI], 1));
        }

        NonOrthogCoeff = max(Foam::asin(sinAlpha)*180.0/M_PI);
    }

    // ZT, 12/Nov/2010
    reduce(NonOrthogCoeff, maxOp<scalar>());

    if (debug)
    {
        Info<< "edgeInterpolation::makeCorrectionVectors() : "
            << "non-orthogonality coefficient = " << NonOrthogCoeff << " deg."
            << endl;
    }

    if (NonOrthogCoeff < 0.1)
    {
        orthogonal_ = true;
        deleteDemandDrivenData(correctionVectors_);
    }
    else
    {
        orthogonal_ = false;
    }

    if (debug)
    {
        Info<< "edgeInterpolation::makeCorrectionVectors() : "
            << "Finished constructing non-orthogonal correction vectors"
            << endl;
    }
}


void edgeInterpolation::makeSkewCorrectionVectors() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeSkewCorrectionVectors() : "
            << "Constructing skew correction vectors"
            << endl;
    }

    skewCorrectionVectors_ = new edgeVectorField
    (
        IOobject
        (
            "skewCorrectionVectors",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimless
    );
    edgeVectorField& SkewCorrVecs = *skewCorrectionVectors_;

    // Set local references to mesh data
    const areaVectorField& C = mesh().areaCentres();
    const edgeVectorField& Ce = mesh().edgeCentres();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const pointField& points = mesh().points();
    const edgeList& edges = mesh().edges();


    forAll(neighbour, edgeI)
    {
        vector P = C[owner[edgeI]];
        vector N = C[neighbour[edgeI]];
        vector S = points[edges[edgeI].start()];
        vector e = edges[edgeI].vec(points);

        scalar alpha = - ( ( (N - P)^(S - P) )&( (N - P)^e ) )/
            ( ( (N - P)^e )&( (N - P)^e ) );

        vector E = S + alpha*e;

        SkewCorrVecs[edgeI] = Ce[edgeI] - E;
    }


    forAll(SkewCorrVecs.boundaryField(), patchI)
    {
        faePatchVectorField& patchSkewCorrVecs =
            SkewCorrVecs.boundaryField()[patchI];

        if (patchSkewCorrVecs.coupled())
        {
            const unallocLabelList& edgeFaces = 
                mesh().boundary()[patchI].edgeFaces();

            const edgeList::subList patchEdges =
                mesh().boundary()[patchI].patchSlice(edges);

            vectorField ngbC = 
                C.boundaryField()[patchI].patchNeighbourField();

            forAll (patchSkewCorrVecs, edgeI)
            {
                vector P = C[edgeFaces[edgeI]];
                vector N = ngbC[edgeI];
                vector S = points[patchEdges[edgeI].start()];
                vector e = patchEdges[edgeI].vec(points);

                scalar alpha = - ( ( (N - P)^(S - P) )&( (N - P)^e ) )/
                    ( ( (N - P)^e )&( (N - P)^e ) );

                vector E = S + alpha*e;

                patchSkewCorrVecs[edgeI] = 
                    Ce.boundaryField()[patchI][edgeI] - E;
            }
        }
        else
        {
            patchSkewCorrVecs = vector::zero;
        }
    }


    scalar skewCoeff = 0.0;

    // Calculating PN arc length
    scalarField lPN(owner.size());

    forAll (owner, edgeI)
    {
        lPN[edgeI] =
            mag
            (
                Ce[edgeI] 
              - SkewCorrVecs[edgeI]
              - C[owner[edgeI]]
            )
          + mag
            (
                C[neighbour[edgeI]]
              - Ce[edgeI] 
              + SkewCorrVecs[edgeI]
            );
    }

    if (lPN.size() > 0)
    {
        skewCoeff = max(mag(SkewCorrVecs.internalField())/mag(lPN));
    }

    if (debug)
    {
        Info<< "edgeInterpolation::makeSkewCorrectionVectors() : "
            << "skew coefficient = " << skewCoeff << endl;
    }

    if (skewCoeff < 0.1)
    {
        skew_ = false;
        deleteDemandDrivenData(skewCorrectionVectors_);
    }
    else
    {
        skew_ = true;
    }

    if (debug)
    {
        Info<< "edgeInterpolation::makeSkewCorrectionVectors() : "
            << "Finished constructing skew correction vectors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
