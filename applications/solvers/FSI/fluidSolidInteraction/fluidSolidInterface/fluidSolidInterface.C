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

\*---------------------------------------------------------------------------*/

#include "fluidSolidInterface.H"
#include "volFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "RectangularMatrix.H"
#include "primitivePatchInterpolation.H"
#include "twoDPointCorrector.H"

#include "tetPointFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "tetPolyPatchInterpolation.H"
#include "tetFemMatrices.H"

#include "fixedValuePointPatchFields.H"
#include "ggiInterpolation.H"
#include "IOmanip.H"
#include "dynamicMotionSolverFvMesh.H"
#include "motionSolver.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "movingWallPressureFvPatchScalarField.H"

// #include "faCFD.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidSolidInterface, 0);
//     defineRunTimeSelectionTable(fluidSolidInterface, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fluidSolidInterface::calcCurrentSolidZonePoints() const
{
    // Find global face zones
    if (currentSolidZonePointsPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcCurrentSolidZonePoints() const"
        )
            << "Current solid zone points alarady exist"
                << abort(FatalError);
    }

    currentSolidZonePointsPtr_ =
        new vectorField(stress().currentFaceZonePoints(solidZoneIndex()));
}

void Foam::fluidSolidInterface::calcCurrentSolidZonePatch() const
{
    // Find global face zones
    if (currentSolidZonePatchPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcCurrentSolidZonePatch() const"
        )
            << "Current solid zone patch alarady exists"
                << abort(FatalError);
    }

    currentSolidZonePatchPtr_ =
        new PrimitivePatch<face, List, const pointField&>
        (
            solidMesh().faceZones()[solidZoneIndex_]().localFaces(),
            currentSolidZonePoints()
        );
}

void Foam::fluidSolidInterface::calcFluidToSolidInterpolator() const
{
    // Find global face zones
    if (fluidToSolidPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcFluidToSolidInterpolator() const"
        )
            << "Fluid to solid interpolator already exists"
                << abort(FatalError);
    }

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    // std::shared_ptr<RBFInterpolation> rbf;
    fluidToSolidPtr_ =
        std::shared_ptr<RBFInterpolation>
        (
            new RBFInterpolation(rbfFunction)
            // new RBFInterpolation(rbfFunctionPtr_)
        );

    vectorField solidZoneFaceCentres =
        currentSolidZonePatch().faceCentres();

    vectorField fluidZoneFaceCentres =
        fluidMesh().faceZones()[fluidZoneIndex_]().faceCentres();

    matrix fluidX(fluidZoneFaceCentres.size(), 3);
    matrix solidX(solidZoneFaceCentres.size(), 3);

    forAll(fluidZoneFaceCentres, faceI)
    {
        fluidX(faceI, 0) = fluidZoneFaceCentres[faceI].x();
        fluidX(faceI, 1) = fluidZoneFaceCentres[faceI].y();
        fluidX(faceI, 2) = fluidZoneFaceCentres[faceI].z();
    }

    forAll(solidZoneFaceCentres, faceI)
    {
        solidX(faceI, 0) = solidZoneFaceCentres[faceI].x();
        solidX(faceI, 1) = solidZoneFaceCentres[faceI].y();
        solidX(faceI, 2) = solidZoneFaceCentres[faceI].z();
    }

    fluidToSolidPtr_->compute(fluidX, solidX);

    Info << "Checking fluid-to-solid interpolator" << endl;
    {
        vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        vectorField fluidZoneFaceCentres
        (
            fluidMesh().faceZones()[fluidZoneIndex_].size(),
            vector::zero
        );

        const label fluidPatchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex_].start();

        forAll (fluidPatchFaceCentres, i)
        {
            fluidZoneFaceCentres
            [
                fluidMesh().faceZones()[fluidZoneIndex_].whichFace
                (
                    fluidPatchStart + i
                )
            ] =
                fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

        vectorField solidPatchFaceCentres =
            vectorField
            (
                solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
            );

        matrix fluidX(fluidZoneFaceCentres.size(), 3);
        // matrix solidX(solidPatchFaceCentres.size(), 3);
        matrix fluidXsolid(solidPatchFaceCentres.size(), 3);

        forAll(fluidZoneFaceCentres, faceI)
        {
            fluidX(faceI, 0) = fluidZoneFaceCentres[faceI].x();
            fluidX(faceI, 1) = fluidZoneFaceCentres[faceI].y();
            fluidX(faceI, 2) = fluidZoneFaceCentres[faceI].z();
        }

        // forAll(solidPatchFaceCentres, faceI)
        // {
        //     solidX(faceI, 0) = solidPatchFaceCentres[faceI].x();
        //     solidX(faceI, 1) = solidPatchFaceCentres[faceI].y();
        //     solidX(faceI, 2) = solidPatchFaceCentres[faceI].z();
        // }

        // fluidToSolidPtr_->compute(fluidX, solidX);
        fluidToSolidPtr_->interpolate(fluidX, fluidXsolid);

        vectorField fluidPatchFaceCentresAtSolid
        (
            solidPatchFaceCentres.size(),
            vector::zero
        );

        forAll(fluidPatchFaceCentresAtSolid, faceI)
        {
            fluidPatchFaceCentresAtSolid[faceI].x() = fluidXsolid(faceI, 0);
            fluidPatchFaceCentresAtSolid[faceI].y() = fluidXsolid(faceI, 1);
            fluidPatchFaceCentresAtSolid[faceI].z() = fluidXsolid(faceI, 2);
        }

        scalar maxDist = gMax
        (
            mag
            (
                 fluidPatchFaceCentresAtSolid
               - solidPatchFaceCentres
            )
        );

        Info << "Fluid-to-solid face interpolation error: " << maxDist
            << endl;
    }
}

// void Foam::fluidSolidInterface::calcGgiFluidToSolidInterpolator() const
// {
//     // Find global face zones
//     if (ggiFluidToSolidPtr_)
//     {
//         FatalErrorIn
//         (
//             "void fluidSolidInterface::"
//             "calcGgiFluidToSolidInterpolator() const"
//         )
//             << "Ggi fluid to solid interpolator already exists"
//                 << abort(FatalError);
//     }

//     ggiFluidToSolidPtr_ =
//         new ggiZoneInterpolation
//         (
//             fluidMesh().faceZones()[fluidZoneIndex_](),
//             solidMesh().faceZones()[solidZoneIndex_](),
//             tensorField(0),
//             tensorField(0),
//             vectorField(0), // Slave-to-master separation. Bug fix
//             0,              // Non-overlapping face tolerances
//             0,              // HJ, 24/Oct/2008
//             true,           // Rescale weighting factors.  Bug fix, MB.
//             ggiInterpolation::AABB //BB_OCTREE  // Octree search, MB.
//         );


//     Info << "Checking fluid-to-solid interpolator" << endl;
//     {
//         vectorField fluidPatchFaceCentres =
//             vectorField
//             (
//                 fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
//             );

//         vectorField fluidZoneFaceCentres
//         (
//             fluidMesh().faceZones()[fluidZoneIndex_].size(),
//             vector::zero
//         );

//         const label fluidPatchStart =
//             fluidMesh().boundaryMesh()[fluidPatchIndex_].start();

//         forAll (fluidPatchFaceCentres, i)
//         {
//             fluidZoneFaceCentres
//             [
//                 fluidMesh().faceZones()[fluidZoneIndex_].whichFace
//                 (
//                     fluidPatchStart + i
//                 )
//             ] =
//                 fluidPatchFaceCentres[i];
//         }

//         // Parallel data exchange: collect faceCentres field on all processors
//         reduce(fluidZoneFaceCentres, sumOp<vectorField>());

//         vectorField solidZoneFaceCentres =
//             ggiFluidToSolidPtr_->masterToSlave
//             (
//                 fluidZoneFaceCentres
//             );

//         vectorField solidPatchFaceCentres
//         (
//             solidMesh().boundaryMesh()[solidPatchIndex_].size(),
//             vector::zero
//         );

//         const label solidPatchStart =
//             solidMesh().boundaryMesh()[solidPatchIndex_].start();

//         forAll(solidPatchFaceCentres, i)
//         {
//             solidPatchFaceCentres[i] =
//                 solidZoneFaceCentres
//                 [
//                     solidMesh().faceZones()[solidZoneIndex_]
//                    .whichFace(solidPatchStart + i)
//                 ];
//         }

//         scalar maxDist = gMax
//         (
//             mag
//             (
//                 solidPatchFaceCentres
//               - solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
//             )
//         );

//         Info << "Fluid-to-solid face interpolation error: " << maxDist
//             << endl;
//     }

//     Info << "Number of uncovered master faces: "
//         << ggiFluidToSolidPtr_->uncoveredMasterFaces().size() << endl;

//     Info << "Number of uncovered slave faces: "
//         << ggiFluidToSolidPtr_->uncoveredSlaveFaces().size() << endl;
// }


void Foam::fluidSolidInterface::calcGgiInterpolator() const
{
    // Create extended ggi interpolation
    if (ggiInterpolatorPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcGgiInterpolator() const"
        )
            << "Ggi interpolator already exists"
                << abort(FatalError);
    }

    // Create copy of solid face zone primitive patch in current configuration

    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(currentSolidZonePointsPtr_);

//     currentSolidZonePatch().movePoints(currentSolidZonePoints());

    Info << "Create extended GGI zone-to-zone interpolator" << endl;

    ggiInterpolatorPtr_ =
        new ggiZoneInterpolation
        (
            fluidMesh().faceZones()[fluidZoneIndex_](),
            currentSolidZonePatch(),
//             solidMesh().faceZones()[solidZoneIndex_](),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            true,           // Patch data is complete on all processors
            SMALL,          // Non-overlapping face tolerances
            SMALL,          // HJ, 24/Oct/2008
            true,           // Rescale weighting factors.  Bug fix, MB.
            ggiInterpolation::BB_OCTREE
            // BB_OCTREE AABB THREE_D_DISTANCE
            // Octree search, MB.
        );

//     currentSolidZonePatch().writeVTK("solidZone");
//     fluidMesh().faceZones()[fluidZoneIndex_]().writeVTK("fluidZone");
//     fluidMesh().boundaryMesh()[fluidPatchIndex_].writeVTK
//     (
//         "fluidPatch",
//         fluidMesh().boundaryMesh()[fluidPatchIndex_],
//         fluidMesh().points()
//     );

    Info << "Checking fluid-to-solid face interpolator" << endl;
    {
        vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        vectorField fluidZoneFaceCentres
        (
            fluidMesh().faceZones()[fluidZoneIndex_].size(),
            vector::zero
        );

        const label fluidPatchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex_].start();

        forAll (fluidPatchFaceCentres, i)
        {
            fluidZoneFaceCentres
            [
                fluidMesh().faceZones()[fluidZoneIndex_].whichFace
                (
                    fluidPatchStart + i
                )
            ] =
                fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZoneFaceCentres =
            ggiInterpolatorPtr_->masterToSlave
            (
                fluidZoneFaceCentres
            );

        vectorField solidPatchFaceCentres
        (
            solidMesh().boundaryMesh()[solidPatchIndex_].size(),
            vector::zero
        );

        const label solidPatchStart =
            solidMesh().boundaryMesh()[solidPatchIndex_].start();

        forAll(solidPatchFaceCentres, i)
        {
            solidPatchFaceCentres[i] =
                solidZoneFaceCentres
                [
                    solidMesh().faceZones()[solidZoneIndex_]
                   .whichFace(solidPatchStart + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                solidPatchFaceCentres
              - solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
            )
        );

        Info << "Fluid-to-solid face interpolation error: " << maxDist
            << endl;
    }

    Info << "Checking solid-to-fluid point interpolator (GGI)" << endl;
    {
        vectorField solidZonePoints_ =
            currentSolidZonePoints();
//             solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        vectorField solidZonePoints =
            ggiInterpolatorPtr_->slaveToMasterPointInterpolate
            (
                solidZonePoints_
            );

        vectorField fluidZonePoints =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        scalar maxDist = gMax
        (
            mag
            (
                fluidZonePoints
              - solidZonePoints
            )
        );

        Info << "Solid-to-fluid point interpolation error (GGI): " << maxDist
            << endl;
    }

    Info << "Number of uncovered master faces: "
        << ggiInterpolatorPtr_->uncoveredMasterFaces().size() << endl;

    Info << "Number of uncovered slave faces: "
        << ggiInterpolatorPtr_ ->uncoveredSlaveFaces().size() << endl;

    ggiInterpolatorPtr_->slavePointDistanceToIntersection();
    ggiInterpolatorPtr_->masterPointDistanceToIntersection();
}


void Foam::fluidSolidInterface::calcSolidToFluidInterpolator() const
{
    // Find global face zones
    if (solidToFluidPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcSolidToFluidInterpolator() const"
        )
            << "Solid to fluid interpolator already exists"
                << abort(FatalError);
    }

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    // std::shared_ptr<RBFInterpolation> rbf;
    solidToFluidPtr_ =
        std::shared_ptr<RBFInterpolation>
        (
            new RBFInterpolation(rbfFunction)
            // new RBFInterpolation(rbfFunctionPtr_)
        );

    vectorField solidZonePoints =
        currentSolidZonePatch().localPoints();

    vectorField fluidZonePoints =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

    matrix fluidX(fluidZonePoints.size(), 3);
    matrix solidX(solidZonePoints.size(), 3);

    forAll(fluidZonePoints, faceI)
    {
        fluidX(faceI, 0) = fluidZonePoints[faceI].x();
        fluidX(faceI, 1) = fluidZonePoints[faceI].y();
        fluidX(faceI, 2) = fluidZonePoints[faceI].z();
    }

    forAll(solidZonePoints, faceI)
    {
        solidX(faceI, 0) = solidZonePoints[faceI].x();
        solidX(faceI, 1) = solidZonePoints[faceI].y();
        solidX(faceI, 2) = solidZonePoints[faceI].z();
    }

    solidToFluidPtr_->compute(solidX, fluidX);

    Info << "Checking solid-to-fluid interpolator" << endl;
    {
        matrix fluidPoints(fluidZonePoints.size(), 3);
        matrix solidPoints(solidZonePoints.size(), 3);
        vectorField fluidZonePointsInterp(fluidZonePoints.size(), vector::zero);


        forAll(solidZonePoints, faceI)
        {
            solidPoints(faceI, 0) = solidZonePoints[faceI].x();
            solidPoints(faceI, 1) = solidZonePoints[faceI].y();
            solidPoints(faceI, 2) = solidZonePoints[faceI].z();
        }

        solidToFluidPtr_->interpolate(solidPoints, fluidPoints);

        forAll(fluidZonePoints, faceI)
        {
            fluidZonePointsInterp[faceI].x() = fluidPoints(faceI, 0);
            fluidZonePointsInterp[faceI].y() = fluidPoints(faceI, 1);
            fluidZonePointsInterp[faceI].z() = fluidPoints(faceI, 2);
        }

        scalar maxDist = gMax
        (
            mag
            (
                fluidZonePointsInterp
              - fluidZonePoints
            )
        );

        Info << "Solid-to-fluid point interpolation error: " << maxDist
            << endl;
    }
}


void Foam::fluidSolidInterface::
calcAccumulatedFluidInterfaceDisplacement() const
{
    // Read accumulated displacement
    if (accumulatedFluidInterfaceDisplacementPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcAccumulatedFluidInterfaceDisplacement() const"
        )
            << "Accumulated displacement field already exists"
                << abort(FatalError);
    }

    // Accumulated fluid interface displacement
    IOobject accumulatedFluidInterfaceDisplacementHeader
    (
        "accumulatedFluidInterfaceDisplacement",
        flow().runTime().timeName(),
        fluidMesh(),
        IOobject::MUST_READ
    );

    if(accumulatedFluidInterfaceDisplacementHeader.headerOk())
    {
        Pout << "Reading accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    flow().runTime().timeName(),
                    fluidMesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        Pout << "Creating accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    flow().runTime().timeName(),
                    fluidMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                vectorField
                (
                    fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
                    vector::zero
                )
            );
    }
}


void Foam::fluidSolidInterface::calcMinEdgeLength() const
{
    // Read accumulated displacement
    if (minEdgeLengthPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcMinEdgeLength() const"
        )
            << "Minimal edge lengths already exist"
                << abort(FatalError);
    }

    minEdgeLengthPtr_ =
        new scalarField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            0
        );
    scalarField& minEdgeLength = *minEdgeLengthPtr_;


    const edgeList& edges =
        fluidMesh().faceZones()[fluidZoneIndex_]().edges();

    const vectorField& points =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

    const labelListList& pointEdges =
        fluidMesh().faceZones()[fluidZoneIndex_]().pointEdges();

    forAll(points, pointI)
    {
        const labelList& curPointEdges = pointEdges[pointI];

        scalar minLength = GREAT;

        forAll(curPointEdges, edgeI)
        {
            const edge& curEdge = edges[curPointEdges[edgeI]];

            scalar Le = curEdge.mag(points);

            if (Le < minLength)
            {
                minLength = Le;
            }
        }

        minEdgeLength[pointI] = minLength;
    }

//     Pout << "Min edge length: " << min(minEdgeLength) << endl;
//     Pout << "gMin edge length: " << gMin(minEdgeLength) << endl;
}

Foam::vectorIOField&
Foam::fluidSolidInterface::accumulatedFluidInterfaceDisplacement()
{
    if (!accumulatedFluidInterfaceDisplacementPtr_)
    {
        calcAccumulatedFluidInterfaceDisplacement();
    }

    return *accumulatedFluidInterfaceDisplacementPtr_;
}


const Foam::scalarField& Foam::fluidSolidInterface::minEdgeLength() const
{
    if (!minEdgeLengthPtr_)
    {
        calcMinEdgeLength();
    }

    return *minEdgeLengthPtr_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::fluidSolidInterface
(
    dynamicFvMesh& fMesh,
    fvMesh& sMesh
)
:
    IOdictionary
    (
        IOobject
        (
            "fsiProperties",
            fMesh.time().constant(),
            fMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    fluidMesh_(fMesh),
    flow_(fluidSolver::New(fluidMesh_)),
    solidMesh_(sMesh),
    stress_(solidSolver::New(solidMesh_)),
    solidPatchIndex_(-1),
    solidZoneIndex_(-1),
    fluidPatchIndex_(-1),
    fluidZoneIndex_(-1),
    currentSolidZonePointsPtr_(NULL),
    currentSolidZonePatchPtr_(NULL),
    fluidToSolidPtr_(NULL),
//     ggiFluidToSolidPtr_(NULL),
    ggiInterpolatorPtr_(NULL),
    solidToFluidPtr_(NULL),
    couplingScheme_(lookup("couplingScheme")),
    relaxationFactor_(readScalar(lookup("relaxationFactor"))),
    aitkenRelaxationFactor_(relaxationFactor_),
    outerCorrTolerance_(readScalar(lookup("outerCorrTolerance"))),
    nOuterCorr_(readInt(lookup("nOuterCorr"))),
    coupled_(lookup("coupled")),
    predictor_(false), //(lookup("predictor")),
    rbfInterpolation_(false), //(lookup("predictor")),
    couplingReuse_(readInt(lookup("couplingReuse"))),
    interfaceDeformationLimit_
    (
        readScalar(lookup("interfaceDeformationLimit"))
    ),
    fluidZonePointsDispl_(),
    fluidZonePointsDisplRef_(),
    fluidZonePointsDisplPrev_(),
    solidZonePointsDispl_(),
    solidZonePointsDisplRef_(),
    interfacePointsDispl_(),
    interfacePointsDisplPrev_(),
    solidZonePressure_(),
    solidZoneTraction_
    (
        IOobject
        (
            "solidZoneTraction",
            fMesh.time().timeName(),
            fMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vectorField()
    ),
    solidZoneTractionPrev_(),
    predictedSolidZoneTraction_
    (
        IOobject
        (
            "predictedSolidZoneTraction",
            fMesh.time().timeName(),
            fMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vectorField()
    ),
    residual_(),
    residualPrev_(),
    maxResidualNorm_(0),
    maxIntDisplNorm_(0),
    outerCorr_(0),
//     closedFluidDomain_(lookup("closedFluidDomain")),
//     refPressure_(0),
//     refPressureIncrement_(0),
//     compressibility_(readScalar(lookup("compressibility"))),
    interpolatorUpdateFrequency_
    (
        readInt(lookup("interpolatorUpdateFrequency"))
    ),
    fluidPatchPointsV_(),
    fluidPatchPointsW_(),
    fluidPatchPointsT_(),
    accumulatedFluidInterfaceDisplacementPtr_(NULL),
    minEdgeLengthPtr_(NULL)
{
    // Solid patch index

    word solidPatchName(lookup("solidPatch"));

    polyPatchID solidPatch(solidPatchName, solidMesh().boundaryMesh());

    if (!solidPatch.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Solid patch name " << solidPatchName << " not found."
                << abort(FatalError);
    }

    solidPatchIndex_ = solidPatch.index();


    // Solid face zone index

    word solidZoneName(lookup("solidZone"));

    faceZoneID solidZone
    (
        solidZoneName,
        solidMesh().faceZones()
    );

    if (!solidZone.active())
    {
        FatalErrorIn("")
            << "Solid face zone name " << solidZoneName
                << " not found.  Please check your face zone definition."
                << abort(FatalError);
    }

    solidZoneIndex_ = solidZone.index();


    // Fluid patch index

    word fluidPatchName(lookup("fluidPatch"));

    polyPatchID fluidPatch(fluidPatchName, fluidMesh().boundaryMesh());

    if (!fluidPatch.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Fluid patch name " << fluidPatchName << " not found."
                << abort(FatalError);
    }

    fluidPatchIndex_ = fluidPatch.index();


    // Fluid face zone index

    word fluidZoneName(lookup("fluidZone"));

    faceZoneID fluidZone
    (
        fluidZoneName,
        fluidMesh().faceZones()
    );

    if (!fluidZone.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Fluid face zone name " << fluidZoneName
                << " not found.  Please check your face zone definition."
                << abort(FatalError);
    }

    fluidZoneIndex_ = fluidZone.index();


    // Check coupling scheme
    if
    (
        (couplingScheme_ == "IQN-ILS")
     || (couplingScheme_ == "Aitken")
     || (couplingScheme_ == "FixedRelaxation")
    )
    {
        Info<< "Selecting coupling scheme " << couplingScheme_ << endl;
    }
    else
    {
        FatalErrorIn
        (
            "fluidSolidInterface::fluidSolidInterface(...)"
        )   << "couplingScheme: " << couplingScheme_
            << " is not a valid choice. "
            << "Options are: IQN-ILS, Aitken, FixedRelaxation"
            << abort(FatalError);
    }

    // Initialize solid zone pressure
    solidZonePressure_ =
        scalarField(solidMesh().faceZones()[solidZoneIndex()].size(), 0.0);

    if (solidZoneTraction_.size() == 0)
    {
        solidZoneTraction_ =
            vectorField
            (
                solidMesh().faceZones()[solidZoneIndex_]().size(),
                vector::zero
            );
    }

    solidZoneTractionPrev_ =
        vectorField
        (
            solidMesh().faceZones()[solidZoneIndex_]().size(),
            vector::zero
        );

    if (predictedSolidZoneTraction_.size() == 0)
    {
        predictedSolidZoneTraction_ =
            vectorField
            (
                solidMesh().faceZones()[solidZoneIndex_]().size(),
                vector::zero
            );
    }

    residual_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    if (found("predictor"))
    {
        predictor_ = Switch(lookup("predictor"));
    }

    if (found("rbfInterpolation"))
    {
        rbfInterpolation_ = Switch(lookup("rbfInterpolation"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::~fluidSolidInterface()
{
    deleteDemandDrivenData(currentSolidZonePointsPtr_);
    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    // deleteDemandDrivenData(fluidToSolidPtr_);
//     deleteDemandDrivenData(ggiFluidToSolidPtr_);
    deleteDemandDrivenData(ggiInterpolatorPtr_);
    // deleteDemandDrivenData(solidToFluidPtr_);
    deleteDemandDrivenData(accumulatedFluidInterfaceDisplacementPtr_);
    deleteDemandDrivenData(minEdgeLengthPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField&
Foam::fluidSolidInterface::currentSolidZonePoints() const
{
    if (!currentSolidZonePointsPtr_)
    {
        calcCurrentSolidZonePoints();
    }

    return *currentSolidZonePointsPtr_;
}

const Foam::PrimitivePatch<Foam::face, Foam::List, const Foam::pointField&>&
Foam::fluidSolidInterface::currentSolidZonePatch() const
{
    if (!currentSolidZonePatchPtr_)
    {
        calcCurrentSolidZonePatch();
    }

    return *currentSolidZonePatchPtr_;
}

const std::shared_ptr<RBFInterpolation>&
Foam::fluidSolidInterface::fluidToSolid() const
{
    if (!fluidToSolidPtr_)
    {
        calcFluidToSolidInterpolator();
    }

    return fluidToSolidPtr_;
}

// const Foam::ggiZoneInterpolation&
// Foam::fluidSolidInterface::ggiFluidToSolid() const
// {
//     if (!ggiFluidToSolidPtr_)
//     {
//         calcGgiFluidToSolidInterpolator();
//     }

//     return *ggiFluidToSolidPtr_;
// }

const Foam::ggiZoneInterpolation&
Foam::fluidSolidInterface::ggiInterpolator() const
{
    if (!ggiInterpolatorPtr_)
    {
        calcGgiInterpolator();
    }

    return *ggiInterpolatorPtr_;
}

const std::shared_ptr<RBFInterpolation>&
Foam::fluidSolidInterface::solidToFluid() const
{
    if (!solidToFluidPtr_)
    {
        calcSolidToFluidInterpolator();
    }

    return solidToFluidPtr_;
}

void Foam::fluidSolidInterface::initializeFields()
{
    fluidZonePointsDispl_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    fluidZonePointsDisplRef_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    fluidZonePointsDisplPrev_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZonePointsDispl_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZonePointsDisplRef_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    residualPrev_ = residual_;
//         vectorField
//         (
//             fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
//             vector::zero
//         );

    residual_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    maxResidualNorm_ = 0;

    outerCorr_ = 0;

    nOuterCorr_ = readInt(lookup("nOuterCorr"));

    outerCorrTolerance_ = readScalar(lookup("outerCorrTolerance"));

    coupled_ = Switch(lookup("coupled"));

    couplingReuse_ = readInt(lookup("couplingReuse"));

    relaxationFactor_ = readScalar(lookup("relaxationFactor"));


    interfacePointsDispl_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    interfacePointsDisplPrev_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

//     refPressure_ += refPressureIncrement_;

//     if (timeIndex_ < runTime().timeIndex())
//     {
//         timeIndex_ = runTime().timeIndex();

//     }
}


void Foam::fluidSolidInterface::updateInterpolator()
{
//     label interpolatorUpdateFrequency_ = 2;

//     Info << runTime().timeIndex() << endl;

    if (!ggiInterpolatorPtr_)
    {
        ggiInterpolator();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex()-1)%interpolatorUpdateFrequency_) == 0)
        {
            deleteDemandDrivenData(ggiInterpolatorPtr_);
            ggiInterpolator();
        }
    }
//     else
//     {
//         if ((runTime().timeIndex()-1) == 0)
//         {
//             deleteDemandDrivenData(ggiInterpolatorPtr_);
//             ggiInterpolator();
//         }
//     }
}


void Foam::fluidSolidInterface::updateDisplacement()
{
    Info << "\nTime = " << flow().runTime().timeName()
        << ", iteration: " << outerCorr() << endl;

//     if (outerCorr_ == 1)
//     {
//         // Cancel residual from previous time step
//         fluidZonePointsDisplPrev() = fluidZonePointsDispl();
//         fluidZonePointsDispl() += residualPrev();
//     }


    if (couplingScheme() == "FixedRelaxation")
    {
        Info << "Current fsi under-relaxation factor: "
            << relaxationFactor() << endl;

        fluidZonePointsDisplPrev() = fluidZonePointsDispl();

        // if (outerCorr_ == 1)
        // {
        //     // Cancel residual from previous time step
        //     fluidZonePointsDispl() += residualPrev();
        // }

        fluidZonePointsDispl() += relaxationFactor()*residual();
    }
    else if (couplingScheme() == "Aitken")
    {
        if (outerCorr_ < 3)
        {
            Info << "Current fsi under-relaxation factor: "
                << relaxationFactor() << endl;

            fluidZonePointsDisplPrev() = fluidZonePointsDispl();

            // if (outerCorr_ == 1)
            // {
            //     // Cancel residual from previous time step
            //     fluidZonePointsDispl() += residualPrev();
            // }

            if ((outerCorr_ == 1) && predictor())
            {
                fluidZonePointsDispl() += residual();
            }
            else
            {
                fluidZonePointsDispl() += relaxationFactor()*residual();
            }
//             fluidZonePointsDispl() += relaxationFactor()*residual();
        }
        else
        {
            aitkenRelaxationFactor() =
               -aitkenRelaxationFactor()
               *(
                    sum
                    (
                        residualPrev()
                      & (residual() - residualPrev())
                    )
                   /(
                        sum
                        (
                            (residual() - residualPrev())
                          & (residual() - residualPrev())
                        )
                    )
                );

            if (Pstream::parRun())
            {
                if(!Pstream::master())
                {
                    aitkenRelaxationFactor() = 0.0;
                }

                //- pass to all procs
                reduce(aitkenRelaxationFactor(), sumOp<scalar>());
            }

            aitkenRelaxationFactor() = mag(aitkenRelaxationFactor());

            if (aitkenRelaxationFactor()>1)
            {
                aitkenRelaxationFactor() = relaxationFactor();
            }

            Info << "Current fsi under-relaxation factor (Aitken): "
                << aitkenRelaxationFactor() << endl;

            fluidZonePointsDisplPrev() = fluidZonePointsDispl();

            fluidZonePointsDispl() +=
                aitkenRelaxationFactor()*residual();
        }
    }
    else if (couplingScheme() == "IQN-ILS")
    {
//      A fluid structure interaction solver with IQN-ILS
//      coupling algorithm (J. Degroote, K.-J. Bathe and J. Vierendeels.
//      Performance of a new partitioned procedure versus a monolithic
//      procedure in fluid-structure interaction. Computers & Structures

        // IQN-ILS
        if (outerCorr_ == 1)
        {
            // Clean up data from old time steps

            Info << "Modes before clean-up : " << fluidPatchPointsT_.size();

            while (true)
            {
                if (fluidPatchPointsT_.size())
                {
                    if
                    (
                        flow().runTime().timeIndex()-couplingReuse()
                      > fluidPatchPointsT_[0]
                    )
                    {
                        for (label i = 0; i < fluidPatchPointsT_.size()-1; i++)
                        {
                            fluidPatchPointsT_[i] = fluidPatchPointsT_[i+1];
                            fluidPatchPointsV_[i] = fluidPatchPointsV_[i+1];
                            fluidPatchPointsW_[i] = fluidPatchPointsW_[i+1];
                        }
                        fluidPatchPointsT_.remove();
                        fluidPatchPointsV_.remove();
                        fluidPatchPointsW_.remove();
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }

            Info << ", modes after clean-up : "
                << fluidPatchPointsT_.size() << endl;
        }
        else if (outerCorr_ == 2)
        {
            // Set reference in the first coupling iteration
            solidZonePointsDisplRef() = solidZonePointsDispl();
            fluidZonePointsDisplRef() = fluidZonePointsDispl();
        }
        else
        {
            // Reference has been set in the first coupling iteration
            fluidPatchPointsV_.append
            (
                (
                    solidZonePointsDispl()
                  - fluidZonePointsDispl()
                )
              - (
                    solidZonePointsDisplRef()
                  - fluidZonePointsDisplRef()
                )
            );

            fluidPatchPointsW_.append
            (
                solidZonePointsDispl()
              - solidZonePointsDisplRef()
            );

            fluidPatchPointsT_.append
            (
                flow().runTime().timeIndex()
            );
        }

        if (fluidPatchPointsT_.size() > 1)
        {
            updateDisplacementUsingIQNILS();
        }
        else
        {
            // Relax the interface displacement
            Info << "Current fsi under-relaxation factor: "
                << relaxationFactor() << endl;

            fluidZonePointsDisplPrev() = fluidZonePointsDispl();

            // if (outerCorr_ == 1)
            // {
            //     // Cancel residual from previous time step
            //     fluidZonePointsDispl() += residualPrev();
            // }

            if ((outerCorr_ == 1) && predictor())
            {
                fluidZonePointsDispl() += residual();
            }
            else
            {
                fluidZonePointsDispl() += relaxationFactor()*residual();
            }
        }
    }


    // Set interface acceleration
    if
    (
        isA<movingWallPressureFvPatchScalarField>
        (
            flow().p().boundaryField()[fluidPatchIndex_]
        )
    )
    {
        Info << "Setting acceleration at fluid side of the interface" << endl;

        vectorField solidZoneAcceleration =
            stress().faceZoneAcceleration
            (
                solidZoneIndex(),
                solidPatchIndex()
            );

        vectorField fluidZoneAcceleration =
            ggiInterpolator().slaveToMaster
            (
                solidZoneAcceleration
            );

        vectorField fluidPatchAcceleration
        (
            fluidMesh().boundary()[fluidPatchIndex()].size(),
            vector::zero
        );

        const label patchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].start();

        forAll(fluidPatchAcceleration, i)
        {
            fluidPatchAcceleration[i] =
                fluidZoneAcceleration
                [
                    fluidMesh().faceZones()[fluidZoneIndex()]
                   .whichFace(patchStart + i)
                ];
        }

        const_cast<movingWallPressureFvPatchScalarField&>
        (
            refCast<const movingWallPressureFvPatchScalarField>
            (
                flow().p().boundaryField()[fluidPatchIndex_]
            )
        ).prevAcceleration() = fluidPatchAcceleration;
    }

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    if (Pstream::parRun())
    {
        if(!Pstream::master())
        {
            fluidZonePointsDispl() *= 0.0;
        }

        //- pass to all procs
        reduce(fluidZonePointsDispl(), sumOp<vectorField>());

        label globalFluidZoneIndex =
            findIndex(flow().globalFaceZones(), fluidZoneIndex());

        if (globalFluidZoneIndex == -1)
        {
            FatalErrorIn
            (
                "fluidSolidInterface::updateDisplacement()"
            ) << "global zone point map is not availabel"
                << abort(FatalError);
        }

        const labelList& map =
            flow().globalToLocalFaceZonePointMap()[globalFluidZoneIndex];

        if(!Pstream::master())
        {
            vectorField fluidZonePointsDisplGlobal =
                fluidZonePointsDispl();

            forAll(fluidZonePointsDisplGlobal, globalPointI)
            {
                label localPoint = map[globalPointI];

                fluidZonePointsDispl()[localPoint] =
                    fluidZonePointsDisplGlobal[globalPointI];
            }
        }
    }
}


Foam::scalar Foam::fluidSolidInterface::updateWeakDisplacement()
{
    vectorField solidZonePointsDisplAtSolid =
        stress().faceZonePointDisplacementIncrement(solidZoneIndex());

    solidZonePointsDispl() =
        ggiInterpolator().slaveToMasterPointInterpolate
        (
            solidZonePointsDisplAtSolid
        );

    vectorField solidZonePointsTotDisplAtSolid =
        stress().faceZonePointDisplacement(solidZoneIndex());

    vectorField solidZonePointsTotDispl =
        ggiInterpolator().slaveToMasterPointInterpolate
        (
            solidZonePointsTotDisplAtSolid
        );


    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

//     Info << "residual: " << average(mag(residual())) << endl;

    fluidZonePointsDisplPrev() = fluidZonePointsDispl();

    fluidZonePointsDispl() += residual();

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    if (Pstream::parRun())
    {
        if(!Pstream::master())
        {
            fluidZonePointsDispl() *= 0.0;
        }

        //- pass to all procs
        reduce(fluidZonePointsDispl(), sumOp<vectorField>());

        label globalFluidZoneIndex =
            findIndex(flow().globalFaceZones(), fluidZoneIndex());

        if (globalFluidZoneIndex == -1)
        {
            FatalErrorIn
            (
                "fluidSolidInterface::updateDisplacement()"
            )   << "global zone point map is not availabel"
                << abort(FatalError);
        }

        const labelList& map =
            flow().globalToLocalFaceZonePointMap()[globalFluidZoneIndex];

        if(!Pstream::master())
        {
            vectorField fluidZonePointsDisplGlobal =
                fluidZonePointsDispl();

            forAll(fluidZonePointsDisplGlobal, globalPointI)
            {
                label localPoint = map[globalPointI];

                fluidZonePointsDispl()[localPoint] =
                    fluidZonePointsDisplGlobal[globalPointI];
            }
        }
    }


    // Set interface acceleration
    if
    (
        isA<elasticWallPressureFvPatchScalarField>
        (
            flow().p().boundaryField()[fluidPatchIndex_]
        )
    )
    {
        vectorField solidZoneAcceleration =
            stress().faceZoneAcceleration
            (
                solidZoneIndex(),
                solidPatchIndex()
            );

        vectorField fluidZoneAcceleration =
            ggiInterpolator().slaveToMaster
            (
                solidZoneAcceleration
            );

        vectorField fluidPatchAcceleration
        (
            fluidMesh().boundary()[fluidPatchIndex()].size(),
            vector::zero
        );

        const label patchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].start();

        forAll(fluidPatchAcceleration, i)
        {
            fluidPatchAcceleration[i] =
                fluidZoneAcceleration
                [
                    fluidMesh().faceZones()[fluidZoneIndex()]
                   .whichFace(patchStart + i)
                ];
        }

        const_cast<elasticWallPressureFvPatchScalarField&>
        (
            refCast<const elasticWallPressureFvPatchScalarField>
            (
                flow().p().boundaryField()[fluidPatchIndex_]
            )
        ).prevAcceleration() = fluidPatchAcceleration;
    }


    // Calculate residual norm

    scalar residualNorm = ::sqrt(gSum(magSqr(residual())));
    scalar residualNorm_2 = residualNorm;

    if (residualNorm > maxResidualNorm_)
    {
        maxResidualNorm_ = residualNorm;
    }

    residualNorm /= maxResidualNorm_ + SMALL;

    Info << "Current fsi relative residual norm: " << residualNorm << endl;

    interfacePointsDisplPrev_ = interfacePointsDispl_;

    interfacePointsDispl_ = solidZonePointsDispl();

    vectorField intTotDispl =
        interfacePointsDispl_ + solidZonePointsTotDispl;
    scalar intTotDisplNorm = ::sqrt(gSum(magSqr(intTotDispl)));
    if (intTotDisplNorm > maxIntDisplNorm_)
    {
        maxIntDisplNorm_ = intTotDisplNorm;
    }

    residualNorm_2 /= maxIntDisplNorm_ + SMALL;

    Info << "Alternative fsi residual: " << residualNorm_2 << endl;

    return min(residualNorm_2, residualNorm);
//     return residualNorm;
}


void Foam::fluidSolidInterface::updateDisplacementUsingIQNILS()
{
    // Consider fluidPatchPointsV as a matrix V
    // with as columns the items
    // in the DynamicList and calculate the QR-decomposition of V
    // with modified Gram-Schmidt
    label cols = fluidPatchPointsV_.size();
    RectangularMatrix<scalar> R(cols, cols, 0.0);
    RectangularMatrix<scalar> C(cols, 1);
    RectangularMatrix<scalar> Rcolsum(1, cols);
    DynamicList<vectorField> Q;

    for (label i = 0; i < cols; i++)
    {
        Q.append(fluidPatchPointsV_[cols-1-i]);
    }

    for (label i = 0; i < cols; i++)
    {
        // Normalize column i
        R[i][i] = Foam::sqrt(sum(Q[i] & Q[i]));
        Q[i] /= R[i][i];

        // Orthogonalize columns to the right of column i
        for (label j = i+1; j < cols; j++)
        {
            R[i][j] = sum(Q[i] & Q[j]);
            Q[j] -= R[i][j]*Q[i];
        }

        // Project minus the residual vector on the Q
        C[i][0] = sum
        (
            Q[i]
          & (
                fluidZonePointsDispl()
              - solidZonePointsDispl()
            )
        );
    }

    // Solve the upper triangular system
    for (label j = 0; j < cols; j++)
    {
        Rcolsum[0][j] = 0.0;
        for (label i = 0; i < j+1; i++)
        {
            Rcolsum[0][j] += cmptMag(R[i][j]);
        }
    }
    scalar epsilon = 1.0E-10*max(Rcolsum);
    for (label i = 0; i < cols; i++)
    {
        if (cmptMag(R[i][i]) > epsilon)
        {
            for (label j = i+1; j < cols; j++)
            {
                R[i][j] /= R[i][i];
            }
            C[i][0] /= R[i][i];
            R[i][i] = 1.0;
        }
    }
    for (label j = cols-1; j >= 0; j--)
    {
        if (cmptMag(R[j][j]) > epsilon)
        {
            for (label i = 0; i < j; i++)
            {
                C[i][0] -= C[j][0]*R[i][j];
            }
        }
        else
        {
            C[j][0] = 0.0;
        }
    }

    fluidZonePointsDisplPrev() = fluidZonePointsDispl();

    fluidZonePointsDispl() = solidZonePointsDispl();

    for (label i = 0; i < cols; i++)
    {
        fluidZonePointsDispl() += fluidPatchPointsW_[i]*C[cols-1-i][0];
    }
}

void Foam::fluidSolidInterface::moveFluidMesh()
{
    // Get fluid patch displacement from fluid zone displacement

    vectorField fluidPatchPointsDispl
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
        vector::zero
    );

    vectorField fluidPatchPointsDisplPrev
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
        vector::zero
    );

    const labelList& fluidPatchMeshPoints =
        fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

    forAll(fluidPatchPointsDispl, pointI)
    {
        label curMeshPointID = fluidPatchMeshPoints[pointI];

        label curFluidZonePointID =
            fluidMesh().faceZones()[fluidZoneIndex()]()
           .whichPoint(curMeshPointID);

        fluidPatchPointsDispl[pointI] =
            fluidZonePointsDispl()[curFluidZonePointID];

        fluidPatchPointsDisplPrev[pointI] =
            fluidZonePointsDisplPrev()[curFluidZonePointID];
    }

    // Move fluid mesh
    const vectorField& n =
        fluidMesh().boundaryMesh()[fluidPatchIndex()].pointNormals();

    primitivePatchInterpolation patchInterpolator
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()]
    );

    scalarField pointDeltaCoeffs =
        patchInterpolator.faceToPointInterpolate
        (
            fluidMesh().boundary()[fluidPatchIndex()].deltaCoeffs()
        );

    scalar delta =
        gMax
        (
            mag
            (
                n
              & (
                    accumulatedFluidInterfaceDisplacement()
                  + fluidPatchPointsDispl
                  - fluidPatchPointsDisplPrev
                )
            )
           *pointDeltaCoeffs
        );

    Info << "Maximal accumulated displacement of interface points: "
        << delta << endl;

    if (false)
//     if (delta < interfaceDeformationLimit())
    {
        // Move only interface points
        pointField newPoints = fluidMesh().allPoints();

        const labelList& meshPoints =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

        forAll (fluidPatchPointsDispl, pointI)
        {
            newPoints[meshPoints[pointI]] +=
                fluidPatchPointsDispl[pointI]
              - fluidPatchPointsDisplPrev[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh_.movePoints(newPoints);

        // Accumulate interface points displacement
        accumulatedFluidInterfaceDisplacement() +=
            fluidPatchPointsDispl
          - fluidPatchPointsDisplPrev;
    }
    else
    {
//         // Move whole fluid mesh
//         const pointField& oldAllPoints = fluidMesh().allPoints();
//         pointField newPoints = fluidMesh().allPoints();

//         const labelList& meshPoints =
//             fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

//         forAll (accumulatedFluidInterfaceDisplacement(), pointI)
//         {
//             newPoints[meshPoints[pointI]] -=
//                 accumulatedFluidInterfaceDisplacement()[pointI];
//         }

//         twoDPointCorrector twoDCorrector(fluidMesh());

//         twoDCorrector.correctPoints(newPoints);

//         fluidMesh_.movePoints(newPoints);

//         accumulatedFluidInterfaceDisplacement() +=
//             fluidPatchPointsDispl
//           - fluidPatchPointsDisplPrev;

        // Check mesh motion solver type
        bool feMotionSolver =
            fluidMesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );

        bool fvMotionSolver =
            fluidMesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    fluidMesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUFluidPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[fluidPatchIndex()]
                );

            tetPolyPatchInterpolation tppi
            (
                refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
            );

            motionUFluidPatch ==
                tppi.pointToPointInterpolate
                (
                    accumulatedFluidInterfaceDisplacement()
                   /flow().runTime().deltaT().value()
                );
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    fluidMesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUFluidPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[fluidPatchIndex()]
                );

            motionUFluidPatch ==
            (
                fluidPatchPointsDispl
              - fluidPatchPointsDisplPrev
            )/flow().runTime().deltaT().value();
        }
        else
        {
            FatalErrorIn("fluidSolidInterface::moveFluidMesh()")
                << "Problem with fluid mesh motion solver selection"
                    << abort(FatalError);
        }

        fluidMesh_.update();

        accumulatedFluidInterfaceDisplacement() =
            vectorField
            (
                accumulatedFluidInterfaceDisplacement().size(),
                vector::zero
            );
    }


    // Move unused fluid mesh points
    {
        vectorField newPoints = fluidMesh().allPoints();

        const labelList& fluidZoneMeshPoints =
            fluidMesh().faceZones()[fluidZoneIndex()]().meshPoints();

        forAll(fluidZonePointsDispl(), pointI)
        {
            if (fluidZoneMeshPoints[pointI] >= fluidMesh().nPoints())
            {
                newPoints[fluidZoneMeshPoints[pointI]] +=
                    fluidZonePointsDispl()[pointI]
                  - fluidZonePointsDisplPrev()[pointI];
            }
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh_.movePoints(newPoints);
    }

//     // Evaluate interface velocity
//     const_cast<volVectorField&>(flow().U())
//    .boundaryField()[fluidPatchIndex()].evaluate();
}


void Foam::fluidSolidInterface::updateForce()
{
    Info << "Setting traction on solid patch" << endl;

    vectorField fluidZoneTraction =
        flow().faceZoneViscousForce
        (
            fluidZoneIndex(),
            fluidPatchIndex()
        );

    // Info << "Viscous force, max: " << max(mag(fluidZoneTraction)) << endl;

    scalarField fluidZonePressure =
        flow().faceZonePressureForce(fluidZoneIndex(), fluidPatchIndex());


    // Fluid zone face normals

    vectorField p =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

    const labelList& mp =
        fluidMesh().faceZones()[fluidZoneIndex_]().meshPoints();
    const vectorField& allPoints = fluidMesh().allPoints();
    forAll(mp, pI)
    {
        p[pI] = allPoints[mp[pI]];
    }

    const faceList& f =
        fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

    vectorField n(f.size(), vector::zero);
    forAll(n, faceI)
    {
        n[faceI] = f[faceI].normal(p);
        n[faceI] /= mag(n[faceI]);
    }

    // Fluid zone total traction
    vectorField fluidZoneTotalTraction =
        fluidZoneTraction - fluidZonePressure*n;

//     vectorField solidZoneTraction =
//         ggiInterpolator().masterToSlave
//         (
//            -fluidZoneTraction
//         );

    vectorField solidZoneTotalTraction
    (
        solidMesh().faceZones()[solidZoneIndex()].size(),
        vector::zero
    );

    if (rbfInterpolation_)
    {
        Info << "... using RBF interpolation" << endl;

        matrix fluidForce(fluidZoneTotalTraction.size(), 3);
        matrix solidForce(solidZoneTotalTraction.size(), 3);

        forAll(fluidZoneTotalTraction, faceI)
        {
            fluidForce(faceI, 0) = fluidZoneTotalTraction[faceI].x();
            fluidForce(faceI, 1) = fluidZoneTotalTraction[faceI].y();
            fluidForce(faceI, 2) = fluidZoneTotalTraction[faceI].z();
        }

        fluidToSolid()->interpolate(fluidForce, solidForce);

        forAll(solidZoneTotalTraction, faceI)
        {
            solidZoneTotalTraction[faceI].x() = -solidForce(faceI, 0);
            solidZoneTotalTraction[faceI].y() = -solidForce(faceI, 1);
            solidZoneTotalTraction[faceI].z() = -solidForce(faceI, 2);
        }
    }
    else
    {
        solidZoneTotalTraction =
            ggiInterpolator().masterToSlave
            (
              - fluidZoneTotalTraction
            );
    }

    solidZonePressure_ =
        ggiInterpolator().masterToSlave
        (
            fluidZonePressure
        );

    if (false)
    {
        volVectorField fluidTraction
        (
            IOobject
            (
                "fluidTraction",
                runTime().timeName(),
                fluidMesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidMesh(),
            dimensionedVector("0", dimPressure, vector::zero)
        );
        fluidTraction.boundaryField()[fluidPatchIndex_] =
          fluidZoneTotalTraction;
        fluidTraction.write();

        volVectorField solidTraction
        (
            IOobject
            (
                "solidTraction",
                runTime().timeName(),
                solidMesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidMesh(),
            dimensionedVector("0", dimPressure, vector::zero)
        );
        solidTraction.boundaryField()[solidPatchIndex_] =
            solidZoneTotalTraction;
        solidTraction.write();
    }

    if (coupled())
    {
//         stress().setPressure
//         (
//             solidPatchIndex(),
//             solidZoneIndex(),
//             solidZonePressure_
//         );

        stress().setTraction
        (
            solidPatchIndex(),
            solidZoneIndex(),
            solidZoneTotalTraction
        );

        // Set interface pressure for elasticWallPressure
        // boundary condition
        if
        (
            isA<elasticWallPressureFvPatchScalarField>
            (
                flow().p().boundaryField()[fluidPatchIndex_]
            )
        )
        {
            const_cast<elasticWallPressureFvPatchScalarField&>
            (
                refCast<const elasticWallPressureFvPatchScalarField>
                (
                    flow().p().boundaryField()[fluidPatchIndex_]
                )
            ).prevPressure() = flow().patchPressureForce(fluidPatchIndex_);
        }
    }
    else
    {
        // Set interface pressure for elasticWallPressure
        // boundary condition
        if
        (
            isA<elasticWallPressureFvPatchScalarField>
            (
                flow().p().boundaryField()[fluidPatchIndex_]
            )
        )
        {
            const_cast<elasticWallPressureFvPatchScalarField&>
            (
                refCast<const elasticWallPressureFvPatchScalarField>
                (
                    flow().p().boundaryField()[fluidPatchIndex_]
                )
            ).prevPressure() = 0;
        }
    }

    // Total force at the fluid side of the interface
    if (true)
    {
        vectorField p =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        const labelList& mp =
            fluidMesh().faceZones()[fluidZoneIndex_]().meshPoints();
        const vectorField& allPoints = fluidMesh().allPoints();
        forAll(mp, pI)
        {
            p[pI] = allPoints[mp[pI]];
        }

        const faceList& f =
            fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        vectorField C(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
            C[faceI] = f[faceI].centre(p);
        }

        vector totalTractionForce =
            sum(fluidZoneTotalTraction*mag(S));

//         vector totalPressureForce = sum(fluidZonePressure*S);

//         Info << setprecision(12);
        Info << "Total force (fluid) = "
            << totalTractionForce << endl;

        // Info << average(C) << ", " << sum(mag(S))
        //     << ", " << sum(fluidZoneTotalTraction) << endl;

//         Info << "Total pressure force (fluid) = "
//             << totalPressureForce << endl;
    }

    // Totla force at the solid side of the interface
    if (true)
    {
        vectorField p =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        const labelList& mp =
            solidMesh().faceZones()[solidZoneIndex_]().meshPoints();
        const vectorField& allPoints = solidMesh().allPoints();
        forAll(mp, pI)
        {
            p[pI] = allPoints[mp[pI]];
        }

        const faceList& f =
            solidMesh().faceZones()[solidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);
        vectorField C(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
            C[faceI] = f[faceI].centre(p);
        }

        vector totalTractionForce =
            sum(solidZoneTotalTraction*mag(S));

//         vector totalPressureForce =
//             sum(solidZonePressure_*S);

//         Info << setprecision(12);
        Info << "Total force (solid) = "
            << totalTractionForce << endl;

        // Info << average(C) << ", " << sum(mag(S))
        //     << ", " << sum(fluidZoneTotalTraction) << endl;

//         Info << sum(mag(S)) << ", " << sum(fluidZoneTotalTraction) << endl;

//         Info << "Total pressure force (solid) = "
//             << totalPressureForce << endl;
    }
}


void Foam::fluidSolidInterface::updateWeakForce()
{
    Info << "Setting weak traction on solid patch" << endl;

    vectorField fluidZoneTraction =
        flow().faceZoneViscousForce
        (
            fluidZoneIndex(),
            fluidPatchIndex()
        );

    scalarField fluidZonePressure =
        flow().faceZonePressureForce
        (
            fluidZoneIndex(),
            fluidPatchIndex()
        );


    // Calculate current interface zone face normals
    // (this is maybe not necessary but I want to be sure
    // that normal is correct for whole zone)
    vectorField n(fluidZonePressure.size(), vector::zero);
    {
        vectorField p =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        const labelList& mp =
            fluidMesh().faceZones()[fluidZoneIndex_]().meshPoints();

        const vectorField& allPoints = fluidMesh().allPoints();

        forAll(mp, pI)
        {
            p[pI] = allPoints[mp[pI]];
        }

        const faceList& f =
            fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

        forAll(n, faceI)
        {
            n[faceI] = f[faceI].normal(p);
            n[faceI] /= mag(n[faceI]);
        }
    }


    // Velocity at fluid and solid side of the interface
    // needed for kinamtically coupled beta scheme

//     // After fluid model solution (Robin bc)
//     vectorField fluidZoneVelocity =
//         flow().faceZoneVelocity
//         (
//             fluidZoneIndex(),
//             fluidPatchIndex()
//         );

//     // Velocity calculated by solving structure model
//     // in previous time step
//     vectorField solidZoneVelocity =
//         stress().faceZoneVelocity
//         (
//             solidZoneIndex(),
//             solidPatchIndex()
//         );

//     vectorField solidZoneVelocityAtFluid =
//         ggiInterpolator().slaveToMaster
//         (
//             solidZoneVelocity
//         );


//     // Pressure force correction according to
//     // kinematically coupled beta scheme
//     scalar rhoSolid = stress().rheology().rho()()[0];
//     scalar mu = stress().rheology().mu()()[0];
//     scalar lambda = stress().rheology().lambda()()[0];
//     scalar ap = sqrt((lambda+2*mu)/rhoSolid);
//     scalar hSolid = ap*flow().runTime().deltaT().value();

//     Info << hSolid << endl;

//     // Fluid interface velocity from previous step
//     vectorField fluidZoneVelocity =
//         flow().U().boundaryField()[fluidPatchIndex_];

//     vectorField fluidZoneVelocityOld =
//         flow().U().oldTime().boundaryField()[fluidPatchIndex_];

//     scalarField fluidZonePressureCorr =
//         rhoSolid*hSolid*
//         ((fluidZoneVelocity - fluidZoneVelocityOld)&n)
//        /flow().runTime().deltaT().value();
//     scalarField fluidZonePressureCorr =
//         rhoSolid*hSolid*
//         ((fluidZoneVelocity - solidZoneVelocityAtFluid)&n)
//        /flow().runTime().deltaT().value();



    // Fluid zone total traction

//     Info << max(fluidZonePressureCorr) << ", "
//         << average(fluidZonePressureCorr) << endl;

    vectorField fluidZoneTotalTraction =
        fluidZoneTraction - fluidZonePressure*n;

    vectorField solidZoneTotalTraction =
        ggiInterpolator().masterToSlave
        (
           -fluidZoneTotalTraction
        );

//     vectorField solidZoneVelocity =
//         ggiInterpolator().masterToSlave
//         (
//             fluidZoneVelocity
//         );

    vectorField solidZoneNormal =
        ggiInterpolator().masterToSlave(-n);

    solidZonePressure_ =
        ggiInterpolator().masterToSlave
        (
            fluidZonePressure
        );

    if (coupled())
    {
//         stress().setVelocityAndTraction
//         (
//             solidPatchIndex(),
//             solidZoneIndex(),
//             solidZoneTotalTraction,
//             solidZoneVelocity,
//             solidZoneNormal
//         );

        stress().setTraction
        (
            solidPatchIndex(),
            solidZoneIndex(),
            solidZoneTotalTraction
        );
    }

    // Total force at the fluid side of the interface
    {
        vectorField p =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        const labelList& mp =
            fluidMesh().faceZones()[fluidZoneIndex_]().meshPoints();
        const vectorField& allPoints = fluidMesh().allPoints();
        forAll(mp, pI)
        {
            p[pI] = allPoints[mp[pI]];
        }

        const faceList& f =
            fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);
//         vectorField C(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
//             C[faceI] = f[faceI].centre(p);
        }

        vector totalTractionForce = sum(fluidZoneTotalTraction*mag(S));

//         Info << setprecision(12);
        Info << "Total force (fluid) = "
            << totalTractionForce << endl;
    }

    // Totla force at the solid side of the interface
    {
        vectorField p =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        const labelList& mp =
            solidMesh().faceZones()[solidZoneIndex_]().meshPoints();
        const vectorField& allPoints = solidMesh().allPoints();
        forAll(mp, pI)
        {
            p[pI] = allPoints[mp[pI]];
        }

        const faceList& f =
            solidMesh().faceZones()[solidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);
//         vectorField C(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
//             C[faceI] = f[faceI].centre(p);
        }

        vector totalTractionForce =
            sum(solidZoneTotalTraction*mag(S));

        //         Info << setprecision(12);
        Info << "Total force (solid) = "
            << totalTractionForce << endl;
    }
}


void Foam::fluidSolidInterface::updateWeakTraction()
{
    Info << "Update weak traction on solid patch" << endl;

    // Calc fluid traction

    const vectorField& p =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();
    const faceList& f =
        fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

    vectorField n(f.size(), vector::zero);
    forAll(n, faceI)
    {
        n[faceI] = f[faceI].normal(p);
        n[faceI] /= mag(n[faceI]);
    }

    vectorField fluidZoneTraction =
        flow().faceZoneViscousForce
        (
            fluidZoneIndex(),
            fluidPatchIndex()
        )
      - flow().faceZonePressureForce(fluidZoneIndex(), fluidPatchIndex())*n;

    vectorField fluidZoneTractionAtSolid =
        ggiInterpolator().masterToSlave
        (
            -fluidZoneTraction
        );

//     scalar beta_ = 1.0;
    scalar beta_ = relaxationFactor_;

    solidZoneTractionPrev_ = solidZoneTraction_;

    solidZoneTraction_ =
        beta_*fluidZoneTractionAtSolid
      + (1.0-beta_)*predictedSolidZoneTraction_;

    predictedSolidZoneTraction_ =
        2*solidZoneTraction_ - solidZoneTractionPrev_;

    if (coupled())
    {
        Info << "Setting weak traction on solid patch" << endl;

        stress().setTraction
        (
            solidPatchIndex(),
            solidZoneIndex(),
            predictedSolidZoneTraction_
        );
    }

    // Total force at the fluid side of the interface
    {
        const vectorField& p =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        const faceList& f =
            fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        vector totalTractionForce = sum(fluidZoneTraction*mag(S));

        Info << "Total force (fluid) = "
            << totalTractionForce << endl;
    }

    // Totla force at the solid side of the interface
    {
        const vectorField& p =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        const faceList& f =
            solidMesh().faceZones()[solidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        vector totalTractionForce =
            sum(fluidZoneTractionAtSolid*mag(S));

        Info << "Total force (solid) = "
            << totalTractionForce << endl;
    }
}


void Foam::fluidSolidInterface::predictAndUpdateForce()
{
    if (coupled())
    {
        Info << "Setting traction on solid patch using prediction" << endl;

        stress().setPressure
        (
            solidPatchIndex(),
            solidZoneIndex(),
            stress().predictPressure(solidPatchIndex(), solidZoneIndex())
        );

        stress().setTraction
        (
            solidPatchIndex(),
            solidZoneIndex(),
            stress().predictTraction(solidPatchIndex(), solidZoneIndex())
        );
    }
}


void Foam::fluidSolidInterface::evolveStress()
{
//     if (closedFluidDomain())
//     {
//         DynamicList<scalar> p0;
//         DynamicList<scalar> dV;

//         scalar requiredVolumeIncrement = 0;
//         forAll(flow().U().boundaryField(), patchI)
//         {
//             if (patchI != fluidPatchIndex())
//             {
//                 requiredVolumeIncrement +=
//                     sum
//                     (
//                         flow().U().boundaryField()[patchI]
//                       & fluidMesh().Sf().boundaryField()[patchI]
//                     )
//                    *runTime().deltaT().value();
//             }
//         }

//         scalar volumeIncrement = 0;

//         label nIter = 0;

//         do
//         {
//             // Calc ref. pressure increment

//             if (dV.size() == 1)
//             {
//                 refPressureIncrement_ += compressibility_*dV[0];
//             }
//             else if (dV.size() >= 2)
//             {
//                 label i = p0.size() - 1;

//                 refPressureIncrement_ =
//                     p0[i-1]
//                   - dV[i-1]*(p0[i] - p0[i-1])
//                    /(dV[i] - dV[i-1] + SMALL);
//             }

//             p0.append(refPressureIncrement_);

//             stress().setPressure
//             (
//                 solidPatchIndex(),
//                 solidZoneIndex(),
//                 solidZonePressure_ + refPressure() + refPressureIncrement()
//             );


//             // Solve solid
//             stress().evolve();

//             // Calculate volume increment

//             const labelList& meshPoints =
//                 solidMesh().faceZones()[solidZoneIndex()]().meshPoints();

//             const faceList& localFaces =
//                 solidMesh().faceZones()[solidZoneIndex()]().localFaces();

//             const vectorField& localPoints =
//                 solidMesh().faceZones()[solidZoneIndex()]().localPoints();


//             vectorField oldLocalPoints =
//                 localPoints
//               + vectorField
//                 (
//                     stress().pointD().oldTime(),
//                     meshPoints
//                 );

//             vectorField newLocalPoints =
//                 localPoints
//               + vectorField
//                 (
//                     stress().pointD(),
//                     meshPoints
//                 );

//             volumeIncrement = 0;

//             forAll(localFaces, faceI)
//             {
//                 volumeIncrement +=
//                     localFaces[faceI].sweptVol
//                     (
//                         oldLocalPoints,
//                         newLocalPoints
//                     );
//             }

//             volumeIncrement -= requiredVolumeIncrement;

//             dV.append(volumeIncrement);
//         }
//         while(mag(volumeIncrement) > SMALL*10 && ++nIter < 10);


//         Info << "Solid volume increment: " << volumeIncrement << endl;
//         Info << "Ref pressure: " << refPressure_ << endl;
//         Info << "Ref pressure increment: " << refPressureIncrement_ << endl;
//         Info << "Calculated compressibility = "
//             << (refPressureIncrement() - p0[0])/(dV[0] + SMALL) << endl;
//     }
//     else
//     {
//         stress().evolve();
//     }

    stress().evolve();
}


Foam::scalar Foam::fluidSolidInterface::updateResidual()
{
    vectorField solidZonePointsDisplAtSolid =
        stress().faceZonePointDisplacementIncrement(solidZoneIndex());

    vectorField solidZonePointsTotDisplAtSolid =
        stress().faceZonePointDisplacement(solidZoneIndex());

    vectorField solidZonePointsTotDispl
    (
        solidZonePointsDispl().size(),
        vector::zero
    );

    if (rbfInterpolation_)
    {
        Info << "Displacement interpolation using RBF interpolation" << endl;

        matrix fluidDispl(solidZonePointsDispl().size(), 3);
        matrix solidDispl(solidZonePointsDisplAtSolid.size(), 3);

        forAll(solidZonePointsDisplAtSolid, pointI)
        {
            solidDispl(pointI, 0) = solidZonePointsDisplAtSolid[pointI].x();
            solidDispl(pointI, 1) = solidZonePointsDisplAtSolid[pointI].y();
            solidDispl(pointI, 2) = solidZonePointsDisplAtSolid[pointI].z();
        }

        solidToFluid()->interpolate(solidDispl, fluidDispl);

        forAll(solidZonePointsDispl(), pointI)
        {
            solidZonePointsDispl()[pointI].x() = fluidDispl(pointI, 0);
            solidZonePointsDispl()[pointI].y() = fluidDispl(pointI, 1);
            solidZonePointsDispl()[pointI].z() = fluidDispl(pointI, 2);
        }

        // Total displacement
        forAll(solidZonePointsTotDisplAtSolid, pointI)
        {
            solidDispl(pointI, 0) = solidZonePointsTotDisplAtSolid[pointI].x();
            solidDispl(pointI, 1) = solidZonePointsTotDisplAtSolid[pointI].y();
            solidDispl(pointI, 2) = solidZonePointsTotDisplAtSolid[pointI].z();
        }

        solidToFluid()->interpolate(solidDispl, fluidDispl);

        forAll(solidZonePointsTotDispl, pointI)
        {
            solidZonePointsTotDispl[pointI].x() = fluidDispl(pointI, 0);
            solidZonePointsTotDispl[pointI].y() = fluidDispl(pointI, 1);
            solidZonePointsTotDispl[pointI].z() = fluidDispl(pointI, 2);
        }
    }
    else
    {
        solidZonePointsDispl() =
            ggiInterpolator().slaveToMasterPointInterpolate
            (
                solidZonePointsDisplAtSolid
            );

        solidZonePointsTotDispl =
            ggiInterpolator().slaveToMasterPointInterpolate
            (
                solidZonePointsTotDisplAtSolid
            );
    }

    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

//     const scalarField& minLe = minEdgeLength();

//     scalar residualNorm = gMax(mag(residual())/(minLe + SMALL));

    scalar residualNorm = ::sqrt(gSum(magSqr(residual())));
    scalar residualNorm_2 = residualNorm;

//     Info << "Current fsi residual norm: " << residualNorm << endl;

    if (residualNorm > maxResidualNorm_)
    {
        maxResidualNorm_ = residualNorm;
    }

//     Info << "Current fsi max residual norm: " << maxResidualNorm_ << endl;

    residualNorm /= maxResidualNorm_ + SMALL;

    Info << "Current fsi relative residual norm: " << residualNorm << endl;

    interfacePointsDisplPrev_ = interfacePointsDispl_;

    interfacePointsDispl_ = solidZonePointsDispl();
//         0.5*(solidZonePointsDispl() + fluidZonePointsDispl());

    vectorField intTotDispl =
        interfacePointsDispl_ + solidZonePointsTotDispl;
    scalar intTotDisplNorm = ::sqrt(gSum(magSqr(intTotDispl)));
    if (intTotDisplNorm > maxIntDisplNorm_)
    {
        maxIntDisplNorm_ = intTotDisplNorm;
    }

    residualNorm_2 /= maxIntDisplNorm_ + SMALL;

//     scalar alterResidual =
//         max(mag(interfacePointsDispl_ - interfacePointsDisplPrev_))
//       /(max(mag(interfacePointsDispl_ + solidZonePointsTotDispl)) + SMALL);

    Info << "Alternative fsi residual: " << residualNorm_2 << endl;

    return min(residualNorm_2, residualNorm);
//     return residualNorm;
}


// bool Foam::fluidSolidInterface::read()
// {
//     if (regIOobject::read())
//     {
//         fluidProperties_ = subDict(type() + "Coeffs");

//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }


// ************************************************************************* //
