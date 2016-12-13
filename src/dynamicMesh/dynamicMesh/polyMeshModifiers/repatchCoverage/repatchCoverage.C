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

#include "repatchCoverage.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "primitiveFacePatch.H"
#include "GGIInterpolationTemplate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(repatchCoverage, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        repatchCoverage,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::repatchCoverage::checkDefinition()
{
    if
    (
        !masterCoveredPatchID_.active()
     || !masterUncoveredPatchID_.active()
     || !slaveCoveredPatchID_.active()
     || !slaveUncoveredPatchID_.active()
     || repatchThreshold_ < 0 || repatchThreshold_ > 1
    )
    {
        FatalErrorIn
        (
            "void Foam::repatchCoverage::checkDefinition()"
        )   << "Not all zones and patches needed in the definition "
            << "have been found.  Please check your mesh definition."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::repatchCoverage::repatchCoverage
(
    const word& name,
    const label index,
    const polyTopoChanger& mme,
    const word& masterCoveredPatchName,
    const word& masterUncoveredPatchName,
    const word& slaveCoveredPatchName,
    const word& slaveUncoveredPatchName,
    const scalar repatchThreshold
)
:
    polyMeshModifier(name, index, mme, true),
    masterCoveredPatchID_(masterCoveredPatchName, mme.mesh().boundaryMesh()),
    masterUncoveredPatchID_
    (
        masterUncoveredPatchName,
        mme.mesh().boundaryMesh()
    ),
    slaveCoveredPatchID_(slaveCoveredPatchName, mme.mesh().boundaryMesh()),
    slaveUncoveredPatchID_
    (
        slaveUncoveredPatchName,
        mme.mesh().boundaryMesh()
    ),
    repatchThreshold_
    (
        Foam::max(SMALL, Foam::min(repatchThreshold, 1 - SMALL))
    ),
    uncMaster_(),
    uncSlave_()
{
    checkDefinition();
}


// Construct from components
Foam::repatchCoverage::repatchCoverage
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    polyMeshModifier(name, index, mme, Switch(dict.lookup("active"))),
    masterCoveredPatchID_
    (
        dict.lookup("masterCoveredPatchName"),
        mme.mesh().boundaryMesh()
    ),
    masterUncoveredPatchID_
    (
        dict.lookup("masterUncoveredPatchName"),
        mme.mesh().boundaryMesh()
    ),
    slaveCoveredPatchID_
    (
        dict.lookup("slaveCoveredPatchName"),
        mme.mesh().boundaryMesh()
    ),
    slaveUncoveredPatchID_
    (
        dict.lookup("slaveUncoveredPatchName"),
        mme.mesh().boundaryMesh()
    ),
    repatchThreshold_
    (
        Foam::max
        (
            SMALL,
            Foam::min(readScalar(dict.lookup("repatchThreshold")), 1 - SMALL)
        )
    ),
    uncMaster_(),
    uncSlave_()
{
    checkDefinition();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::repatchCoverage::~repatchCoverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::repatchCoverage::changeTopology() const
{
    // Check that masks are empty
    if (!uncMaster_.empty() || !uncSlave_.empty())
    {
        FatalErrorIn("bool repatchCoverage::changeTopology() const")
            << "Uncovered masks are not empty.  Topo change is out of sync"
            << abort(FatalError);
    }

    // Get mesh reference
    const polyMesh& mesh = topoChanger().mesh();

    const pointField& allPoints = mesh.allPoints();
    const faceList& allFaces = mesh.allFaces();

    // Collect all faces for interpolation

    // Master side
    const polyPatch& masterCoveredPatch =
        mesh.boundaryMesh()[masterCoveredPatchID_.index()];

    const polyPatch& masterUncoveredPatch =
        mesh.boundaryMesh()[masterUncoveredPatchID_.index()];

    faceList masterSide
    (
        masterCoveredPatch.size()
      + masterUncoveredPatch.size()
    );

    {
        label nFaces = 0;

        // Insert covered faces
        for
        (
            label faceI = masterCoveredPatch.start();
            faceI < masterCoveredPatch.start() + masterCoveredPatch.size();
            faceI++
        )
        {
            masterSide[nFaces] = allFaces[faceI];
            nFaces++;
        }

        // Insert uncovered faces
        for
        (
            label faceI = masterUncoveredPatch.start();
            faceI < masterUncoveredPatch.start() + masterUncoveredPatch.size();
            faceI++
        )
        {
            masterSide[nFaces] = allFaces[faceI];
            nFaces++;
        }
    }

    // Slave side
    const polyPatch& slaveCoveredPatch =
        mesh.boundaryMesh()[slaveCoveredPatchID_.index()];

    const polyPatch& slaveUncoveredPatch =
        mesh.boundaryMesh()[slaveUncoveredPatchID_.index()];

    faceList slaveSide
    (
        slaveCoveredPatch.size()
      + slaveUncoveredPatch.size()
    );

    {
        label nFaces = 0;

        // Insert covered faces
        for
        (
            label faceI = slaveCoveredPatch.start();
            faceI < slaveCoveredPatch.start() + slaveCoveredPatch.size();
            faceI++
        )
        {
            slaveSide[nFaces] = allFaces[faceI];
            nFaces++;
        }

        // Insert uncovered faces
        for
        (
            label faceI = slaveUncoveredPatch.start();
            faceI < slaveUncoveredPatch.start() + slaveUncoveredPatch.size();
            faceI++
        )
        {
            slaveSide[nFaces] = allFaces[faceI];
            nFaces++;
        }
    }

    // Create interpolator
    primitiveFacePatch master(masterSide, allPoints);
    primitiveFacePatch slave(slaveSide, allPoints);

    if
    (
        Pstream::parRun()
     && (
            (masterSide.empty() && !slaveSide.empty())
         || (!masterSide.empty() && slaveSide.empty())
        )
    )
    {
        FatalErrorIn("bool repatchCoverage::changeTopology() const")
            << "Parallel run with partial visibility"
            << abort(FatalError);
    }

    if (masterSide.empty() && slaveSide.empty())
    {
        // Empty master and slave patches.  No local topo change
        return false;
    }

    GGIInterpolation<primitiveFacePatch, primitiveFacePatch> patchToPatch
    (
        master,
        slave,
        tensorField::zero,  // forwardT
        tensorField::zero,  // reverseT
        vectorField::zero,  // separation
        false               // Patches are not global
    );

    // Check uncovered master and slave faces

    // Check uncovered master faces
    // All faces between 0 and masterCoveredPatch.size() - 1 should be covered
    // All faces between masterCoveredPatch.size() - 1 and
    //      uncoveredMaster.size() - 1 should be uncovered
    // If this is not the case, faces should be changed

    // Master side
    label nMasterChanges = 0;

    {
        // Set size of uncovered master mask.  This indicates a topo change
        // is in preparation
        uncMaster_.setSize(master.size(), false);

        const labelList& umf = patchToPatch.uncoveredMasterFaces();

        forAll(umf, umfI)
        {
            uncMaster_[umf[umfI]] = true;
        }

        // Check coverage
        forAll (uncMaster_, mfI)
        {
            if (uncMaster_[mfI] && mfI < masterCoveredPatch.size())
            {
                // Found uncovered master in covered section
                nMasterChanges++;
            }

            if (!uncMaster_[mfI] && mfI >= masterCoveredPatch.size())
            {
                // Found covered master in uncovered section
                nMasterChanges++;
            }
        }
    }

    // Slave side
    label nSlaveChanges = 0;

    {
        // Set size of uncovered master mask.  This indicates a topo change
        // is in preparation
        uncSlave_.setSize(slave.size(), false);

        const labelList& umf = patchToPatch.uncoveredSlaveFaces();

        forAll(umf, umfI)
        {
            uncSlave_[umf[umfI]] = true;
        }

        // Check coverage
        forAll (uncSlave_, mfI)
        {
            if (uncSlave_[mfI] && mfI < slaveCoveredPatch.size())
            {
                // Found uncovered slave in covered section
                nSlaveChanges++;
            }

            if (!uncSlave_[mfI] && mfI >= slaveCoveredPatch.size())
            {
                // Found covered slave in uncovered section
                nSlaveChanges++;
            }
        }
    }

    if (nMasterChanges > 0 || nSlaveChanges > 0)
    {
        InfoIn("bool repatchCoverage::changeTopology() const")
            << "Changing " << nMasterChanges << " master and "
            << nSlaveChanges << " slave faces"
            << endl;

        return true;
    }
    else
    {
        // Clear uncovered masks to indicate that the topological change
        // is not required
        uncMaster_.clear();
        uncSlave_.clear();

        return false;
    }
}


void Foam::repatchCoverage::setRefinement(polyTopoChange& ref) const
{
    // Get mesh reference
    const polyMesh& mesh = topoChanger().mesh();

    const faceList& allFaces = mesh.allFaces();
    const faceZoneMesh& faceZones = mesh.faceZones();

    const labelList& owner = mesh.faceOwner();

    // Master side
    {
        // Get master patches

        const polyPatch& masterCoveredPatch =
            mesh.boundaryMesh()[masterCoveredPatchID_.index()];

        const label cStart = masterCoveredPatch.start();
        const label cSize = masterCoveredPatch.size();

        const polyPatch& masterUncoveredPatch =
            mesh.boundaryMesh()[masterUncoveredPatchID_.index()];

        const label uncStart = masterUncoveredPatch.start();

        // Adjust coverage
        forAll (uncMaster_, faceI)
        {
            if (uncMaster_[faceI] && faceI < masterCoveredPatch.size())
            {
                // Found uncovered face in covered section

                // Get face index
                const label faceIndex = cStart + faceI;

                // Find zone index
                const label faceZoneIndex = faceZones.whichZone(faceIndex);

                // Face flip
                bool faceZoneFlip = false;

                if (faceZoneIndex > -1)
                {
                    // Find face in zone
                    const label fizIndex =
                        faceZones[faceZoneIndex].whichFace(faceIndex);

                    faceZoneFlip =
                        faceZones[faceZoneIndex].flipMap()[fizIndex];
                }

                // Found uncovered master in covered section
                // Move to uncovered patch
                ref.setAction
                (
                    polyModifyFace
                    (
                        allFaces[faceIndex],             // modified face
                        faceIndex,                       // modified face index
                        owner[faceIndex],                // owner
                        -1,                              // neighbour
                        false,                           // face flip
                        masterUncoveredPatch.index(),    // patch for face
                        false,                           // remove from zone
                        faceZoneIndex,                   // zone for face
                        faceZoneFlip                     // face flip in zone
                    )
                );
            }

            if (!uncMaster_[faceI] && faceI >= masterCoveredPatch.size())
            {
                // Found covered face in uncovered section

                // Get face index
                const label faceIndex = uncStart - cSize + faceI;

                // Find zone index
                const label faceZoneIndex = faceZones.whichZone(faceIndex);

                // Face flip
                bool faceZoneFlip = false;

                if (faceZoneIndex > -1)
                {
                    // Find face in zone
                    const label fizIndex =
                        faceZones[faceZoneIndex].whichFace(faceIndex);

                    faceZoneFlip =
                        faceZones[faceZoneIndex].flipMap()[fizIndex];
                }

                // Found uncovered master in covered section
                // Move to uncovered patch
                ref.setAction
                (
                    polyModifyFace
                    (
                        allFaces[faceIndex],             // modified face
                        faceIndex,                       // modified face index
                        owner[faceIndex],                // owner
                        -1,                              // neighbour
                        false,                           // face flip
                        masterCoveredPatch.index(),      // patch for face
                        false,                           // remove from zone
                        faceZoneIndex,                   // zone for face
                        faceZoneFlip                     // face flip in zone
                    )
                );
            }
        }
    }

    {
        // Get slave patches

        const polyPatch& slaveCoveredPatch =
            mesh.boundaryMesh()[slaveCoveredPatchID_.index()];

        const label cStart = slaveCoveredPatch.start();
        const label cSize = slaveCoveredPatch.size();

        const polyPatch& slaveUncoveredPatch =
            mesh.boundaryMesh()[slaveUncoveredPatchID_.index()];

        const label uncStart = slaveUncoveredPatch.start();

        // Check coverage
        forAll (uncSlave_, faceI)
        {
            if (uncSlave_[faceI] && faceI < slaveCoveredPatch.size())
            {
                // Found uncovered face in covered section

                // Get face index
                const label faceIndex = cStart + faceI;

                // Find zone index
                const label faceZoneIndex = faceZones.whichZone(faceIndex);

                // Face flip
                bool faceZoneFlip = false;

                if (faceZoneIndex > -1)
                {
                    // Find face in zone
                    const label fizIndex =
                        faceZones[faceZoneIndex].whichFace(faceIndex);

                    faceZoneFlip =
                        faceZones[faceZoneIndex].flipMap()[fizIndex];
                }

                // Found uncovered slave in covered section
                // Move to uncovered patch
                ref.setAction
                (
                    polyModifyFace
                    (
                        allFaces[faceIndex],             // modified face
                        faceIndex,                       // modified face index
                        owner[faceIndex],                // owner
                        -1,                              // neighbour
                        false,                           // face flip
                        slaveUncoveredPatch.index(),    // patch for face
                        false,                           // remove from zone
                        faceZoneIndex,                   // zone for face
                        faceZoneFlip                     // face flip in zone
                    )
                );
            }

            if (!uncSlave_[faceI] && faceI >= slaveCoveredPatch.size())
            {
                // Found covered face in uncovered section

                // Get face index
                const label faceIndex = uncStart - cSize + faceI;

                // Find zone index
                const label faceZoneIndex = faceZones.whichZone(faceIndex);

                // Face flip
                bool faceZoneFlip = false;

                if (faceZoneIndex > -1)
                {
                    // Find face in zone
                    const label fizIndex =
                        faceZones[faceZoneIndex].whichFace(faceIndex);

                    faceZoneFlip =
                        faceZones[faceZoneIndex].flipMap()[fizIndex];
                }

                // Found uncovered slave in covered section
                // Move to uncovered patch
                ref.setAction
                (
                    polyModifyFace
                    (
                        allFaces[faceIndex],             // modified face
                        faceIndex,                       // modified face index
                        owner[faceIndex],                // owner
                        -1,                              // neighbour
                        false,                           // face flip
                        slaveCoveredPatch.index(),      // patch for face
                        false,                           // remove from zone
                        faceZoneIndex,                   // zone for face
                        faceZoneFlip                     // face flip in zone
                    )
                );
            }
        }
    }


    // Clear uncovered masks to indicate that the topological change
    // has been performed
    uncMaster_.clear();
    uncSlave_.clear();
}


void Foam::repatchCoverage::modifyMotionPoints(pointField& motionPoints) const
{}


void Foam::repatchCoverage::updateMesh(const mapPolyMesh&)
{
    // Mesh has changed topologically.  Update local topological data
    const polyMesh& mesh = topoChanger().mesh();

    masterCoveredPatchID_.update(mesh.boundaryMesh());
    masterUncoveredPatchID_.update(mesh.boundaryMesh());
    slaveCoveredPatchID_.update(mesh.boundaryMesh());
    slaveUncoveredPatchID_.update(mesh.boundaryMesh());
}


void Foam::repatchCoverage::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name() << nl
        << masterCoveredPatchID_.name() << nl
        << masterUncoveredPatchID_.name() << nl
        << slaveCoveredPatchID_.name() << nl
        << slaveUncoveredPatchID_.name() << endl;
}


void Foam::repatchCoverage::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type()
        << token::END_STATEMENT << nl
        << "    masterCoveredPatchName " << masterCoveredPatchID_.name()
        << token::END_STATEMENT << nl
        << "    masterUncoveredPatchName " << masterUncoveredPatchID_.name()
        << token::END_STATEMENT << nl
        << "    slaveCoveredPatchName " << slaveCoveredPatchID_.name()
        << token::END_STATEMENT << nl
        << "    slaveUncoveredPatchName " << slaveUncoveredPatchID_.name()
        << token::END_STATEMENT << nl
        << "    repatchThreshold " << repatchThreshold_
        << token::END_STATEMENT << nl
        << "    active " << active()
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
