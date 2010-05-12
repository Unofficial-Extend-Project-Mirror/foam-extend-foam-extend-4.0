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

#include "vtkPV3Foam.H"

// Foam includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "vtkPV3FoamReader.h"

// local headers
#include "vtkPV3FoamAddToSelection.H"
#include "vtkPV3FoamUpdateInformationFields.H"

// VTK includes
#include "vtkDataArraySelection.h"


// * * * * * * * * * * * * * * * Private Classes * * * * * * * * * * * * * * //

namespace Foam
{

//- A class for reading zone information without requiring a mesh
class zonesEntries
:
    public regIOobject,
    public PtrList<entry>
{

public:

    // Constructors

        explicit zonesEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream("regIOobject"))
        {
            close();
        }

   // Member functions

        bool writeData(Ostream&) const
        {
            notImplemented("zonesEntries::writeData(Ostream&) const");
            return false;
        }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::vtkPV3Foam::readZoneNames(const word& zoneType)
{
    wordList zoneNames;

    // mesh not loaded - read from file
    IOobject ioObj
    (
        zoneType,
        dbPtr_().findInstance
        (
            polyMesh::meshSubDir,
            zoneType,
            IOobject::READ_IF_PRESENT
        ),
        polyMesh::meshSubDir,
        dbPtr_(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (ioObj.headerOk())
    {
        zonesEntries zones(ioObj);

        zoneNames.setSize(zones.size());
        forAll (zones, zoneI)
        {
            zoneNames[zoneI] = zones[zoneI].keyword();
        }
    }

    return zoneNames;
}


void Foam::vtkPV3Foam::updateInformationInternalMesh()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInformationInternalMesh" << endl;
    }

    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    // Determine number of meshes available
    HashTable<const fvMesh*> meshObjects = dbPtr_().lookupClass<const fvMesh>();
    nMesh_ = meshObjects.size();

    // Determine regions (internal mesh and patches...)
    //- Add internal mesh as first entry
    selectInfoVolume_ = arraySelection->GetNumberOfArrays();
    arraySelection->AddArray("internalMesh");
    selectInfoVolume_ += 1;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInformationInternalMesh" << endl;
    }

}


void Foam::vtkPV3Foam::updateInformationLagrangian()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInformationLagrangian" << nl
            << "    " << dbPtr_->timePath()/"lagrangian" << endl;
    }

    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    // Search for list of lagrangian objects for this time
    fileNameList cloudDirs
    (
        readDir(dbPtr_->timePath()/"lagrangian", fileName::DIRECTORY)
    );

    selectInfoLagrangian_ = arraySelection->GetNumberOfArrays();

    if (cloudDirs.size())
    {
        arraySelection->AddArray("lagrangian");
        selectInfoLagrangian_ += 1;

        Info<< "... added cloudDirs\n";

        if (cloudDirs.size() > 1)
        {
            WarningIn("void Foam::vtkPV3Foam::updateInformationLagrangian()")
                << "Multiple lagrangian clouds identified. Currently only able "
                << "to process ONE cloud: " << cloudDirs[0]
                << endl;
        }
        // Set cloud name to first cloud found
        // TODO - multiple clouds
        cloudName_ = cloudDirs[0];
    }
    else
    {
        if (debug)
        {
            Info<< "... no clouds identified in " <<nl
                << "    " <<dbPtr_->timePath()/"lagrangian" << endl;
        }
    }

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInformationLagrangian" << endl;
    }
}


void Foam::vtkPV3Foam::updateInformationPatches()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInformationPatches"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }

    vtkDataArraySelection *arraySelection = reader_->GetRegionSelection();
    selectInfoPatches_ = arraySelection->GetNumberOfArrays();

    int nPatches = 0;

    if (meshPtr_)
    {
        const polyBoundaryMesh& patches = meshPtr_->boundaryMesh();
        forAll (patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.size())
            {
                // Add patch to GUI region list
                arraySelection->AddArray
                (
                    (pp.name() + " - patch").c_str()
                );

                ++nPatches;
            }
        }
    }
    else
    {
        // Read patches
        polyBoundaryMeshEntries patchEntries
        (
            IOobject
            (
                "boundary",
                dbPtr_().findInstance(polyMesh::meshSubDir, "boundary"),
                polyMesh::meshSubDir,
                dbPtr_(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // Start regions at patches
        forAll (patchEntries, entryI)
        {
            label nFaces
            (
                readLabel(patchEntries[entryI].dict().lookup("nFaces"))
            );

            // Valid patch if nFace > 0
            if (nFaces)
            {
                // Add patch to GUI region list
                arraySelection->AddArray
                (
                    (patchEntries[entryI].keyword() + " - patch").c_str()
                );

                ++nPatches;
            }
        }
    }

    selectInfoPatches_ += nPatches;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInformationPatches" << endl;
    }
}


void Foam::vtkPV3Foam::updateInformationZones()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInformationZones"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }

    vtkDataArraySelection *arraySelection = reader_->GetRegionSelection();

    wordList namesLst;

    //
    // cellZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->cellZones().names();
    }
    else
    {
        namesLst = readZoneNames("cellZones");
    }

    selectInfoCellZones_ = arraySelection->GetNumberOfArrays();
    forAll (namesLst, elemI)
    {
        arraySelection->AddArray((namesLst[elemI] + " - cellZone").c_str());
    }
    selectInfoCellZones_ += namesLst.size();
    zoneSuperCells_.setSize(selectInfoCellZones_.size());


    //
    // faceZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->faceZones().names();
    }
    else
    {
        namesLst = readZoneNames("faceZones");
    }

    selectInfoFaceZones_ = arraySelection->GetNumberOfArrays();
    forAll (namesLst, elemI)
    {
        arraySelection->AddArray((namesLst[elemI] + " - faceZone").c_str());
    }
    selectInfoFaceZones_ += namesLst.size();


    //
    // pointZones information
    // ~~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->pointZones().names();
    }
    else
    {
        namesLst = readZoneNames("pointZones");
    }

    selectInfoPointZones_ = arraySelection->GetNumberOfArrays();
    forAll (namesLst, elemI)
    {
        arraySelection->AddArray((namesLst[elemI] + " - pointZone").c_str());
    }
    selectInfoPointZones_ += namesLst.size();


    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInformationZones" << endl;
    }
}


void Foam::vtkPV3Foam::updateInformationSets()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInformationSets" << endl;
    }

    vtkDataArraySelection *arraySelection = reader_->GetRegionSelection();

    // Add names of sets
    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().findInstance(polyMesh::meshSubDir, "faces"),
        polyMesh::meshSubDir/"sets"
    );


    selectInfoCellSets_ = arraySelection->GetNumberOfArrays();
    selectInfoCellSets_ += addToSelection<cellSet>
    (
        arraySelection,
        objects,
        " - cellSet"
    );
    csetSuperCells_.setSize(selectInfoCellSets_.size());

    selectInfoFaceSets_ = arraySelection->GetNumberOfArrays();
    selectInfoFaceSets_ += addToSelection<faceSet>
    (
        arraySelection,
        objects,
        " - faceSet"
    );

    selectInfoPointSets_ = arraySelection->GetNumberOfArrays();
    selectInfoPointSets_ += addToSelection<pointSet>
    (
        arraySelection,
        objects,
        " - pointSet"
    );

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInformationSets" << endl;
    }
}


void Foam::vtkPV3Foam::updateInformationLagrangianFields()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInformationLagrangianFields"
            << endl;
    }

    vtkDataArraySelection *arraySelection =
        reader_->GetLagrangianFieldSelection();

    // preserve the currently selected values
    const stringList selectedEntries = getSelectedArrayEntries
    (
        arraySelection
    );
    arraySelection->RemoveAllArrays();

    // TODO - currently hard-coded to ONE cloud
    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().timeName(),
        "lagrangian"/cloudName_
    );

    addToSelection<IOField<label> >
    (
        arraySelection,
        objects
    );
    addToSelection<IOField<scalar> >
    (
        arraySelection,
        objects
    );
    addToSelection<IOField<vector> >
    (
        arraySelection,
        objects
    );
    addToSelection<IOField<sphericalTensor> >
    (
        arraySelection,
        objects
    );
    addToSelection<IOField<symmTensor> >
    (
        arraySelection,
        objects
    );
    addToSelection<IOField<tensor> >
    (
        arraySelection,
        objects
    );

    // restore the currently enabled values
    setSelectedArrayEntries
    (
        arraySelection,
        selectedEntries
    );

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateInformationLagrangianFields - "
            << "lagrangian objects.size() = " << objects.size() << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
