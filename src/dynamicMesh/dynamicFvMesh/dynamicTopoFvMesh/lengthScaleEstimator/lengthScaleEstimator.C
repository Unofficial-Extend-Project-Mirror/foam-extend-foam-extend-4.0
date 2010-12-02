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

Class
    lengthScaleEstimator

Description
    Implementation of lengthScaleEstimator members

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "volFields.H"
#include "lengthScaleEstimator.H"
#include "processorPolyPatch.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(lengthScaleEstimator,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyMesh and dictionary
lengthScaleEstimator::lengthScaleEstimator
(
    const polyMesh& mesh
)
:
    mesh_(mesh),
    ratioMin_(0.0),
    ratioMax_(0.0),
    growthFactor_(1.0),
    minLengthScale_(VSMALL),
    maxLengthScale_(GREAT),
    curvatureDeviation_(0.0),
    spatialRes_(-1),
    proximityBins_(0),
    sliceThreshold_(VSMALL),
    sliceHoldOff_(0),
    sliceBoxes_(0),
    field_("none"),
    fieldLength_(0.0),
    lowerRefineLevel_(0.001),
    upperRefineLevel_(0.999),
    meanScale_(-1.0),
    maxRefineLevel_(labelMax)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lengthScaleEstimator::~lengthScaleEstimator()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Check for legitimacy of patches
void lengthScaleEstimator::checkPatches(const wordList& patchList) const
{
    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    forAll(patchList, wordI)
    {
        if (boundary.findPatchID(patchList[wordI]) < 0)
        {
            FatalErrorIn
            (
                "void lengthScaleEstimator::checkPatches"
                "(const wordList& patchList) const"
            )
                << " Could not find patch: "
                << patchList[wordI] << nl
                << abort(FatalError);
        }
    }
}


// Prepare for proximity-based refinement, if necessary
void lengthScaleEstimator::prepareProximityPatches()
{
    if (!proximityPatches_.size())
    {
        return;
    }

    if (debug)
    {
        Info << "Preparing patches for proximity-based refinement...";
    }

    // Clear out existing lists.
    proximityBins_.clear();
    proximityBins_.setSize(997, labelList(0));

    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    bool setSpatialRes = false;

    // Loop through all proximity patches and spatially hash patch faces.
    forAll(boundary, patchI)
    {
        if
        (
            (proximityPatches_.found(boundary[patchI].name())) ||
            (boundary[patchI].type() == "symmetryPlane")
        )
        {
            const polyPatch& proxPatch = boundary[patchI];

            // Construct a bounding-box of face centres.
            // Do not synchronize in parallel, since the patch
            // may not be present on all sub-domains.
            boundBox box(proxPatch.faceCentres(), false);

            const point& bMin = box.min();
            const point& bMax = box.max();

            // Further hashing requires this information.
            proxBoundBox_.min() = bMin;
            proxBoundBox_.max() = bMax;

            // Build a list of face indices
            labelList faceIndices
            (
                identity(proxPatch.size()) + proxPatch.start()
            );

            // For spatial resolution, pick a face on this patch.
            if (!setSpatialRes)
            {
                scalar eLength =
                (
                    Foam::sqrt(mag(proxPatch.faceAreas()[0]))
                );

                spatialRes_ =
                (
                    label(::floor(mag(bMax - bMin) / (3.0 * eLength)))
                );

                setSpatialRes = true;
            }

            spatialHash
            (
                proxPatch.faceCentres(),
                faceIndices,
                box,
                spatialRes_,
                proximityBins_
            );
        }
    }

    if (debug)
    {
        Info << "Done." << endl;
    }
}


// Perform spatial hashing on a set of points
void lengthScaleEstimator::spatialHash
(
    const pointField& pointLocations,
    const labelList& pointIndices,
    const boundBox& box,
    const label resolution,
    labelListList& bins
)
{
    label binSize = bins.size(), nD = resolution;

    const point& bMin = box.min();
    const point& bMax = box.max();

    // Extend bounding-box dimensions a bit to avoid edge-effects.
    scalar ext = 0.02*(mag(bMax - bMin));

    // Define an inverse grid-cell size.
    scalar xL = nD/(bMax.x() - bMin.x() + ext);
    scalar yL = nD/(bMax.y() - bMin.y() + ext);
    scalar zL = nD/(bMax.z() - bMin.z() + ext);

    // Loop through all points and bin them.
    forAll(pointLocations, pointI)
    {
        // Translate to boundBox minimum.
        point p = pointLocations[pointI] - bMin;

        // Hash the position.
        label i = label(mag(::floor(p.x()*xL)));
        label j = label(mag(::floor(p.y()*yL)));
        label k = label(mag(::floor(p.z()*zL)));

        label pos = mag(((k*nD*nD)+(j*nD)+i) % binSize);

        // Store the index.
        label newSize = bins[pos].size() + 1;
        bins[pos].setSize(newSize, pointIndices[pointI]);
    }
}


// Send length-scale info across processors
void lengthScaleEstimator::writeLengthScaleInfo
(
    const labelList& cellLevels,
    const scalarList& lengthScale
)
{
    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    // Clear existing buffers
    sendLblBuffer_.clear();
    recvLblBuffer_.clear();
    sendLvlBuffer_.clear();
    recvLvlBuffer_.clear();
    sendSclBuffer_.clear();
    recvSclBuffer_.clear();

    // Number of sent faces across processors
    labelList nSendFaces(boundary.size(), 0);
    labelList nRecvFaces(boundary.size(), 0);

    // Corresponding face labels
    sendLblBuffer_.setSize(boundary.size());
    sendLvlBuffer_.setSize(boundary.size());
    recvLblBuffer_.setSize(boundary.size());
    recvLvlBuffer_.setSize(boundary.size());

    // Length-scales corresponding to face-labels
    sendSclBuffer_.setSize(boundary.size());
    recvSclBuffer_.setSize(boundary.size());

    // Fill send buffers with cell-level and length-scale info.
    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            const labelList& fCells = boundary[pI].faceCells();

            // Set the initial buffer size
            sendLblBuffer_[pI].setSize(boundary[pI].size(), 0);
            sendLvlBuffer_[pI].setSize(boundary[pI].size(), 0);
            sendSclBuffer_[pI].setSize(boundary[pI].size(), 0.0);

            forAll(fCells, faceI)
            {
                label cI = fCells[faceI];

                // Does the adjacent cell have a non-zero level?
                if (cellLevels[cI] > 0)
                {
                    // Fill the send buffer.
                    sendLblBuffer_[pI][nSendFaces[pI]] = faceI;
                    sendLvlBuffer_[pI][nSendFaces[pI]] = cellLevels[cI];
                    sendSclBuffer_[pI][nSendFaces[pI]] = lengthScale[cI];

                    nSendFaces[pI]++;
                }
            }

            // Resize to actual value
            sendLblBuffer_[pI].setSize(nSendFaces[pI]);
            sendLvlBuffer_[pI].setSize(nSendFaces[pI]);
            sendSclBuffer_[pI].setSize(nSendFaces[pI]);

            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label neiProcNo = pp.neighbProcNo();

            // First perform a blocking send/receive of the number of faces.
            OPstream::write
            (
                Pstream::blocking,
                neiProcNo,
                reinterpret_cast<const char*>(&(nSendFaces[pI])),
                sizeof(label)
            );

            IPstream::read
            (
                Pstream::blocking,
                neiProcNo,
                reinterpret_cast<char*>(&(nRecvFaces[pI])),
                sizeof(label)
            );
        }
    }

    // Send info to neighbouring processors.
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            label neiProcNo = pp.neighbProcNo();

            if (debug > 4)
            {
                Pout << " Processor patch " << patchI << ' ' << pp.name()
                     << " communicating with " << neiProcNo
                     << "  Sending: " << nSendFaces[patchI]
                     << "  Recving: " << nRecvFaces[patchI]
                     << endl;
            }

            // Next, perform a non-blocking send/receive of buffers.
            // But only if the buffer size is non-zero.
            if (nSendFaces[patchI] != 0)
            {
                OPstream::write
                (
                    Pstream::nonBlocking,
                    neiProcNo,
                    reinterpret_cast<const char*>(&(sendLblBuffer_[patchI][0])),
                    sendLblBuffer_[patchI].size()*sizeof(label)
                );

                OPstream::write
                (
                    Pstream::nonBlocking,
                    neiProcNo,
                    reinterpret_cast<const char*>(&(sendLvlBuffer_[patchI][0])),
                    sendLblBuffer_[patchI].size()*sizeof(label)
                );

                OPstream::write
                (
                    Pstream::nonBlocking,
                    neiProcNo,
                    reinterpret_cast<const char*>(&(sendSclBuffer_[patchI][0])),
                    sendSclBuffer_[patchI].size()*sizeof(scalar)
                );
            }

            if (nRecvFaces[patchI] != 0)
            {
                // Size the receive buffers
                recvLblBuffer_[patchI].setSize(nRecvFaces[patchI], 0);
                recvLvlBuffer_[patchI].setSize(nRecvFaces[patchI], 0);
                recvSclBuffer_[patchI].setSize(nRecvFaces[patchI], 0.0);

                IPstream::read
                (
                    Pstream::nonBlocking,
                    neiProcNo,
                    reinterpret_cast<char*>(&(recvLblBuffer_[patchI][0])),
                    nRecvFaces[patchI]*sizeof(label)
                );

                IPstream::read
                (
                    Pstream::nonBlocking,
                    neiProcNo,
                    reinterpret_cast<char*>(&(recvLvlBuffer_[patchI][0])),
                    nRecvFaces[patchI]*sizeof(label)
                );

                IPstream::read
                (
                    Pstream::nonBlocking,
                    neiProcNo,
                    reinterpret_cast<char*>(&(recvSclBuffer_[patchI][0])),
                    nRecvFaces[patchI]*sizeof(scalar)
                );
            }
        }
    }
}


// Receive length-scale info across processors
void lengthScaleEstimator::readLengthScaleInfo
(
    const label level,
    label& visitedCells,
    labelList& cellLevels,
    UList<scalar>& lengthScale,
    labelHashSet& levelCells
) const
{
    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Wait for all transfers to complete.
    OPstream::waitRequests();
    IPstream::waitRequests();

    // Now re-visit cells and update length-scales.
    forAll(boundary, patchI)
    {
        if (recvLblBuffer_[patchI].size())
        {
            const labelList& fCells = boundary[patchI].faceCells();

            forAll(recvLblBuffer_[patchI], i)
            {
                label nLabel = recvLblBuffer_[patchI][i];
                label pLevel = recvLvlBuffer_[patchI][i];

                label cI = fCells[nLabel];

                label& ngbLevel = cellLevels[cI];

                // For an unvisited cell, update the level
                bool unvisited = false;

                if (ngbLevel == 0)
                {
                    ngbLevel = pLevel + 1;
                    unvisited = true;
                }

                // Visit all neighbours and re-calculate length-scale
                const cell& cellCheck = mesh_.cells()[cI];

                scalar sumLength = 0.0;
                label nTouchedNgb = 0;

                forAll(cellCheck, faceI)
                {
                    label sLevel = -1, sLCell = -1;
                    label pF = boundary.whichPatch(cellCheck[faceI]);

                    if (pF == -1)
                    {
                        // Internal face. Determine neighbour.
                        if (own[cellCheck[faceI]] == cI)
                        {
                            sLCell = nei[cellCheck[faceI]];
                        }
                        else
                        {
                            sLCell = own[cellCheck[faceI]];
                        }

                        sLevel = cellLevels[sLCell];

                        if ((sLevel < ngbLevel) && (sLevel > 0))
                        {
                            sumLength += lengthScale[sLCell];

                            nTouchedNgb++;
                        }
                    }
                    else
                    if (isA<processorPolyPatch>(boundary[pF]))
                    {
                        // Determine the local index.
                        label local =
                        (
                            boundary[pF].whichFace(cellCheck[faceI])
                        );

                        // Is this label present in the list?
                        label j = -1;

                        if ((j = findIndex(recvLblBuffer_[pF], local)) > -1)
                        {
                            if
                            (
                                (recvLvlBuffer_[pF][j] < ngbLevel) &&
                                (recvLvlBuffer_[pF][j] > 0)
                            )
                            {
                                sumLength += recvSclBuffer_[pF][j];

                                nTouchedNgb++;
                            }
                        }
                    }
                    else
                    if (!isFreePatch(pF))
                    {
                        sumLength +=
                        (
                            fixedLengthScale(cellCheck[faceI], pF, true)
                        );

                        nTouchedNgb++;
                    }
                }

                sumLength /= nTouchedNgb;

                // Scale the length and assign to this cell
                if (level < maxRefineLevel_)
                {
                    sumLength *= growthFactor_;
                }
                else
                if (meanScale_ > 0.0)
                {
                    // If a mean scale has been specified,
                    // override the value
                    sumLength = meanScale_;
                }

                lengthScale[cI] = sumLength;

                if (unvisited)
                {
                    levelCells.insert(cI);

                    visitedCells++;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Read edge refinement options from the dictionary
void lengthScaleEstimator::readRefinementOptions
(
    const dictionary& refineDict,
    bool reRead,
    bool mandatory
)
{
    // Save old values
    scalar oldRatioMax = ratioMax_;
    scalar oldRatioMin = ratioMin_;
    scalar oldGrowthFactor = growthFactor_;
    scalar oldMinScale = minLengthScale_;
    scalar oldMaxScale = maxLengthScale_;

    ratioMax_ = readScalar(refineDict.lookup("bisectionRatio"));
    ratioMin_ = readScalar(refineDict.lookup("collapseRatio"));
    growthFactor_ = readScalar(refineDict.lookup("growthFactor"));

    // Sanity check: Are ratios and growth-factor correctly specified?
    if
    (
        ratioMin_ > ratioMax_ ||
        ratioMin_ < VSMALL ||
        ratioMax_ < VSMALL ||
        growthFactor_ < 1.0
    )
    {
        FatalErrorIn
        (
            "void lengthScaleEstimator::readRefinementOptions"
            "(const dictionary&, bool, bool)"
        )
            << " Options are incorrectly specified." << nl
            << " ratioMin: " << ratioMin_ << nl
            << " ratioMax: " << ratioMax_ << nl
            << " growthFactor: " << growthFactor_ << nl
            << abort(FatalError);
    }

    if (refineDict.found("minLengthScale") || mandatory)
    {
        minLengthScale_ = readScalar(refineDict.lookup("minLengthScale"));
    }

    if (refineDict.found("maxLengthScale") || mandatory)
    {
        maxLengthScale_ = readScalar(refineDict.lookup("maxLengthScale"));
    }

    // Sanity check: Are min/max length scales correctly specified?
    if (minLengthScale_ > maxLengthScale_)
    {
        FatalErrorIn
        (
            "void lengthScaleEstimator::readRefinementOptions"
            "(const dictionary&, bool, bool)"
        )
            << " Length-scales are incorrectly specified." << nl
            << " minLengthScale: " << minLengthScale_ << nl
            << " maxLengthScale: " << maxLengthScale_ << nl
            << abort(FatalError);
    }

    if (reRead)
    {
        // Check if values have changed, and report it.
        if (mag(oldRatioMax - ratioMax_) > SMALL)
        {
            Info << "\tOld ratioMax: " << oldRatioMax << nl
                 << "\tNew ratioMax: " << ratioMax_ << endl;
        }

        if (mag(oldRatioMin - ratioMin_) > SMALL)
        {
            Info << "\tOld ratioMin: " << oldRatioMin << nl
                 << "\tNew ratioMin: " << ratioMin_ << endl;
        }

        if (mag(oldGrowthFactor - growthFactor_) > SMALL)
        {
            Info << "\tOld growthFactor: " << oldGrowthFactor << nl
                 << "\tNew growthFactor: " << growthFactor_ << endl;
        }

        if (mag(oldMinScale - minLengthScale_) > SMALL)
        {
            Info << "\tOld minLengthScale: " << oldMinScale << nl
                 << "\tNew minLengthScale: " << minLengthScale_ << endl;
        }

        if (mag(oldMaxScale - maxLengthScale_) > SMALL)
        {
            Info << "\tOld maxLengthScale: " << oldMaxScale << nl
                 << "\tNew maxLengthScale: " << maxLengthScale_ << endl;
        }
    }

    if (refineDict.found("fixedLengthScalePatches") || mandatory)
    {
        fixedPatches_ = refineDict.subDict("fixedLengthScalePatches");

        // Ensure that patches are legitimate.
        checkPatches(fixedPatches_.toc());
    }

    if (refineDict.found("freeLengthScalePatches") || mandatory)
    {
        freePatches_ = refineDict.subDict("freeLengthScalePatches");

        // Ensure that patches are legitimate.
        checkPatches(freePatches_.toc());

        // Check if fixed and free patches are conflicting
        if (fixedPatches_.size() && freePatches_.size())
        {
            wordList fixedPatchList = fixedPatches_.toc();
            wordList freePatchList = freePatches_.toc();

            forAll(fixedPatchList, wordI)
            {
                if (findIndex(freePatchList, fixedPatchList[wordI]) > -1)
                {
                    FatalErrorIn
                    (
                        "void lengthScaleEstimator::readRefinementOptions"
                        "(const dictionary&, bool, bool)"
                    )
                        << " Conflicting fixed/free patches." << nl
                        << " Patch: " << fixedPatchList[wordI] << nl
                        << abort(FatalError);
                }
            }
        }
    }

    if (refineDict.found("curvaturePatches") || mandatory)
    {
        curvaturePatches_ = refineDict.subDict("curvaturePatches");

        // Ensure that patches are legitimate.
        checkPatches(curvaturePatches_.toc());

        curvatureDeviation_ =
        (
            readScalar(refineDict.lookup("curvatureDeviation"))
        );

        if (curvatureDeviation_ > 1.0 || curvatureDeviation_ < 0.0)
        {
            FatalErrorIn
            (
                "void lengthScaleEstimator::readRefinementOptions"
                "(const dictionary&, bool, bool)"
            )
                << " Curvature deviation out of range [0..1]"
                << abort(FatalError);
        }
    }

    if (refineDict.found("proximityPatches") || mandatory)
    {
        proximityPatches_ =
        (
            refineDict.subDict("proximityPatches")
        );

        // Ensure that patches are legitimate.
        checkPatches(proximityPatches_.toc());

        // Check if a threshold for slicing has been specified.
        if (refineDict.found("sliceThreshold") || mandatory)
        {
            sliceThreshold_ =
            (
                readScalar(refineDict.lookup("sliceThreshold"))
            );

            // Cap the threshold value
            // sliceThreshold_ = Foam::max(sliceThreshold_, minLengthScale_);
        }
    }

    // Check if edge-refinement is to be avoided on any patches
    if (refineDict.found("noModificationPatches") || mandatory)
    {
        wordList noModPatches =
        (
            refineDict.subDict("noModificationPatches").toc()
        );

        // Ensure that patches are legitimate.
        checkPatches(noModPatches);

        noModPatchIDs_.setSize(noModPatches.size());

        label indexI = 0;

        forAll(noModPatches, wordI)
        {
            const word& patchName = noModPatches[wordI];

            noModPatchIDs_[indexI++] =
            (
                mesh_.boundaryMesh().findPatchID(patchName)
            );
        }
    }

    // Check if refinement is to be performed based on a field-value
    if (refineDict.found("fieldRefinement") || mandatory)
    {
        field_ = word(refineDict.lookup("fieldRefinement"));

        // Lookup a specified length-scale
        fieldLength_ = readScalar(refineDict.lookup("fieldLengthScale"));

        lowerRefineLevel_ = readScalar(refineDict.lookup("lowerRefineLevel"));
        upperRefineLevel_ = readScalar(refineDict.lookup("upperRefineLevel"));

        // Sanity check: Are refinement levels correctly specified?
        if (lowerRefineLevel_ > upperRefineLevel_)
        {
            FatalErrorIn
            (
                "void lengthScaleEstimator::readRefinementOptions"
                "(const dictionary&, bool, bool)"
            )
                << " Refinement levels are incorrectly specified." << nl
                << " lowerRefineLevel: " << lowerRefineLevel_ << nl
                << " upperRefineLevel: " << upperRefineLevel_ << nl
                << abort(FatalError);
        }
    }

    // Check if a max refinement level has been specified
    if (refineDict.found("maxRefineLevel") || mandatory)
    {
        maxRefineLevel_ = readLabel(refineDict.lookup("maxRefineLevel"));
    }

    // Update switch for mean-scale
    if (refineDict.found("meanScale") || mandatory)
    {
        meanScale_ = readScalar(refineDict.lookup("meanScale"));
    }
}


//- Set explicitly coupled patch information
void lengthScaleEstimator::setCoupledPatches
(
    const dictionary& coupledPatches
)
{
    label indexI = 0;

    // Determine master and slave patches
    wordList masterPatches(coupledPatches.size());
    wordList slavePatches(coupledPatches.size());

    forAllConstIter(dictionary, coupledPatches, dIter)
    {
        const dictionary& dictI = dIter().dict();

        masterPatches[indexI] = word(dictI.lookup("master"));
        slavePatches[indexI] = word(dictI.lookup("slave"));

        indexI++;
    }

    // Ensure that patches are legitimate.
    checkPatches(masterPatches);

    // Check whether coupled patches are fixedPatches as well.
    forAll(masterPatches, wordI)
    {
        word pName(masterPatches[wordI]);

        if (fixedPatches_.found(pName))
        {
            // Add the slave patch to the list as well.
            // If it already exists, over-ride the value.
            fixedPatches_.add
            (
                slavePatches[wordI],
                fixedPatches_[pName][0].scalarToken(),
                true
            );
        }
    }
}


//- Calculate the length scale field
void lengthScaleEstimator::calculateLengthScale
(
    UList<scalar>& lengthScale
)
{
    label level = 1, visitedCells = 0;
    labelList cellLevels(mesh_.nCells(), 0);

    // Check for allocation
    if (lengthScale.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "void lengthScaleEstimator::calculateLengthScale"
            "(UList<scalar>& lengthScale)"
        )
            << " Field is incorrectly sized." << nl
            << " Field size: " << lengthScale.size()
            << " nCells: " << mesh_.nCells()
            << abort(FatalError);
    }

    // HashSet to keep track of cells in each level
    labelHashSet levelCells;

    // Prepare for proximity-based refinement, if necessary
    prepareProximityPatches();

    // Obtain the cellCells addressing list
    const labelList& own = mesh_.faceOwner();
    const labelListList& cc = mesh_.cellCells();
    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    forAll(boundary,patchI)
    {
        // Skip floating length-scale patches
        if (isFreePatch(patchI))
        {
            continue;
        }

        const polyPatch& bdyPatch = boundary[patchI];

        label pStart = bdyPatch.start();

        forAll(bdyPatch,faceI)
        {
            label ownCell = own[pStart+faceI];

            if (cellLevels[ownCell] != 0)
            {
                continue;
            }

            cellLevels[ownCell] = level;

            lengthScale[ownCell] =
            (
                fixedLengthScale(pStart+faceI, patchI, true) * growthFactor_
            );

            levelCells.insert(ownCell);

            visitedCells++;
        }
    }

    // If a field has been specified, use that.
    if (field_ != "none")
    {
        const volScalarField& vFld =
        (
            mesh_.objectRegistry::lookupObject<volScalarField>(field_)
        );

        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();

        forAll(nei, faceI)
        {
            label ownCell = own[faceI], neiCell = nei[faceI];

            scalar fAvg = 0.5 * (vFld[ownCell] + vFld[neiCell]);

            if ((fAvg > lowerRefineLevel_) && (fAvg < upperRefineLevel_))
            {
                // Set the scale for cells on either side
                if (!cellLevels[ownCell])
                {
                    cellLevels[ownCell] = level;
                    lengthScale[ownCell] = fieldLength_;

                    levelCells.insert(ownCell);

                    visitedCells++;
                }

                if (!cellLevels[neiCell])
                {
                    cellLevels[neiCell] = level;
                    lengthScale[neiCell] = fieldLength_;

                    levelCells.insert(neiCell);

                    visitedCells++;
                }
            }
        }
    }

    bool doneWithSweeps = false;

    // Perform multiple sweeps through the mesh...
    while (!doneWithSweeps)
    {
        if (Pstream::parRun())
        {
            writeLengthScaleInfo(cellLevels, lengthScale);
        }

        // Loop through cells of the current level
        labelList currLvlCells = levelCells.toc();
        levelCells.clear();

        // Loop through cells, and increment neighbour
        // cells of the current level
        forAll(currLvlCells,cellI)
        {
            // Obtain the cells neighbouring this one
            const labelList& cList = cc[currLvlCells[cellI]];

            forAll(cList, indexI)
            {
                label& ngbLevel = cellLevels[cList[indexI]];

                if (ngbLevel == 0)
                {
                    ngbLevel = level + 1;

                    // Compute the mean of the existing
                    // neighbour length-scales
                    const labelList& ncList = cc[cList[indexI]];
                    scalar sumLength = 0.0;
                    label nTouchedNgb = 0;

                    forAll(ncList, indexJ)
                    {
                        label sLevel = cellLevels[ncList[indexJ]];

                        if ((sLevel < ngbLevel) && (sLevel > 0))
                        {
                            sumLength += lengthScale[ncList[indexJ]];

                            nTouchedNgb++;
                        }
                    }

                    sumLength /= nTouchedNgb;

                    // Scale the length and assign to this cell
                    if (level < maxRefineLevel_)
                    {
                        sumLength *= growthFactor_;
                    }
                    else
                    if (meanScale_ > 0.0)
                    {
                        // If a mean scale has been specified,
                        // override the value
                        sumLength = meanScale_;
                    }

                    lengthScale[cList[indexI]] = sumLength;

                    levelCells.insert(cList[indexI]);

                    visitedCells++;
                }
            }
        }

        if (Pstream::parRun())
        {
            readLengthScaleInfo
            (
                level,
                visitedCells,
                cellLevels,
                lengthScale,
                levelCells
            );
        }

        if (debug > 2)
        {
            Pout << "Processed level: " << level << nl
                 << " Visited: " << visitedCells
                 << " out of " << mesh_.nCells() << endl;
        }

        // Move on to the next level
        level++;

        if (visitedCells >= mesh_.nCells())
        {
            doneWithSweeps = true;
        }

        // Wait for everyone to complete.
        reduce(doneWithSweeps, andOp<bool>());
    }

    if (debug)
    {
        Info << "Max Length Scale: " << maxLengthScale_ << endl;
        Info << "Length Scale sweeps: " << level << endl;
    }

    // Check if everything went okay
    if (visitedCells != mesh_.nCells())
    {
        FatalErrorIn
        (
            "void lengthScaleEstimator::calculateLengthScale"
            "(UList<scalar>& lengthScale)"
        )
            << " Algorithm did not visit every cell in the mesh."
            << " Something's messed up." << nl
            << " Visited cells: " << visitedCells
            << " nCells: " << mesh_.nCells()
            << abort(FatalError);
    }

    // Wait for transfers before continuing.
    if (Pstream::parRun())
    {
        OPstream::waitRequests();
        IPstream::waitRequests();
    }
}


} // End namespace Foam

// ************************************************************************* //
