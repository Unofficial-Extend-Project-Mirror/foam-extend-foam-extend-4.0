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

#include "multiTopoBodyFvMesh.H"
#include "foamTime.H"
#include "mapPolyMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiTopoBodyFvMesh, 0);
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        multiTopoBodyFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiTopoBodyFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if (!topoChanger_.empty())
    {
        Info<< "void multiTopoBodyFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh.  " << bodies_.size()
        << " bodies found" << endl;

    // Guess sizes of lists
    DynamicList<pointZone*> pz(bodies_.size());
    DynamicList<faceZone*> fz(4*bodies_.size() + faceZones().size());
    DynamicList<cellZone*> cz(bodies_.size());

    // Copy the face zones associated with the GGI interfaces
    if (faceZones().size() > 0)
    {
        // Copy point zones
        Info << "Copying existing point zones" << endl;

        forAll (pointZones(), i)
        {
            pz.append(pointZones()[i].clone(pointZones()).ptr());
        }

        forAll (faceZones(), i)
        {
            fz.append(faceZones()[i].clone(faceZones()).ptr());
        }

        forAll (cellZones(), i)
        {
            cz.append(cellZones()[i].clone(cellZones()).ptr());
        }
    }

    Info << "Adding point, face and cell zones" << endl;
    forAll (bodies_, bodyI)
    {
        bodies_[bodyI].addZones(pz, fz, cz);
    }

    {
        List<pointZone*> pzList;
        pzList.transfer(pz.shrink());

        List<faceZone*> fzList;
        fzList.transfer(fz.shrink());

        List<cellZone*> czList;
        czList.transfer(cz.shrink());

        removeZones();
        addZones(pzList, fzList, czList);
    }

    // Add modifiers
    {
        topoChanger_.setSize(4*bodies_.size());
        label nextI = 0;

        forAll (bodies_, bodyI)
        {
            bodies_[bodyI].addModifiers(topoChanger_, nextI);
        }

        Info<< "Adding topology modifiers.  nModifiers = " << nextI << endl;

        // Resize and set topo changer
        topoChanger_.setSize(nextI);
        topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
        topoChanger_.write();
    }

    // Write mesh with new zones, updating GGI zoning etc
    syncUpdateMesh();

    write();
}


Foam::tmp<Foam::vectorField> Foam::multiTopoBodyFvMesh::pointMotion() const
{
    // Accumulate point motion
    tmp<vectorField> tpointMotion
    (
        new vectorField(allPoints().size(), vector::zero)
    );

    vectorField& pointMotion = tpointMotion();

    forAll (bodies_, bodyI)
    {
        pointMotion += bodies_[bodyI].pointMotion();
    }

    return tpointMotion;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::multiTopoBodyFvMesh::multiTopoBodyFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    bodies_()
{
    // Read body entries from the dictionary
    PtrList<entry> bodyEntries(dict_.lookup("bodies"));
    bodies_.setSize(bodyEntries.size());

    forAll (bodyEntries, bodyI)
    {
        bodies_.set
        (
            bodyI,
            new topoBody
            (
                bodyEntries[bodyI].keyword(),
                *this,
                bodyEntries[bodyI].dict()
            )
        );
    }

    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiTopoBodyFvMesh::~multiTopoBodyFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiTopoBodyFvMesh::update()
{
    // Save old points
    pointField oldPointsNew = allPoints();
    pointField newPoints = allPoints() + pointMotion();

    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();
    bool localMeshChanged = topoChangeMap->morphing();
    bool globalMeshChanged = localMeshChanged;
    reduce(globalMeshChanged, orOp<bool>());

    if (globalMeshChanged)
    {
        forAll (bodies_, bodyI)
        {
            bodies_[bodyI].updateTopology();
        }
    }

    if (localMeshChanged)
    {
        // Map old points onto the new mesh
        // Not needed: only a single motion in time-step so old points
        // are already in the correct postion for mesh motion.
        // HJ, 19/May/2014
//         pointField mappedOldPointsNew(allPoints().size());
//         mappedOldPointsNew.map(oldPointsNew, topoChangeMap->pointMap());

//         // Note: using setOldPoints instead of movePoints.
//         // HJ, 23/Aug/2015
//         setOldPoints(mappedOldPointsNew);
//         resetMotion();
//         setV0();

        // Get new points from preMotion
        newPoints = topoChangeMap().preMotionPoints() + pointMotion();
    }
    else
    {
        // Topo change on other processors
    }

    Info << "Executing mesh motion" << endl;
    movePoints(newPoints);

    return localMeshChanged;
}


// ************************************************************************* //
