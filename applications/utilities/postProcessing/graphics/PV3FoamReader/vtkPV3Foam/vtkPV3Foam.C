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
#include "Time.H"
#include "fvMesh.H"
#include "IOobjectList.H"
#include "patchZones.H"
#include "vtkPV3FoamReader.h"
#include "IFstream.H"

// VTK includes
#include "vtkCharArray.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkInformation.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::vtkPV3Foam, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "vtkPV3FoamAddToSelection.H"
#include "vtkPV3FoamUpdateInformationFields.H"


void Foam::vtkPV3Foam::AddToBlock
(
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo,
    vtkDataSet* dataset,
    const string& blockName
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);
    if (blockDO && !block)
    {
        FatalErrorIn("Foam::vtkPV3Foam::AddToBlock")
            << "Block already has a vtkDataSet assigned to it" << nl << endl;
        return;
    }

    if (!block)
    {
        block = vtkMultiBlockDataSet::New();
        output->SetBlock(blockNo, block);
        block->Delete();
    }

    if (block)
    {
        if (debug)
        {
            Info<< "block[" << blockNo << "] has "
                << block->GetNumberOfBlocks()
                <<  " datasets prior to adding set " << datasetNo
                <<  " with name: " << blockName << endl;
        }

        // when assigning dataset 0, also name the parent block
        if (!datasetNo && selector.name())
        {
            output->GetMetaData(blockNo)->Set
            (
                vtkCompositeDataSet::NAME(),
                selector.name()
            );
        }
    }


    block->SetBlock(datasetNo, dataset);

    if (blockName.size())
    {
        block->GetMetaData(datasetNo)->Set
        (
            vtkCompositeDataSet::NAME(), blockName.c_str()
        );
    }

}


vtkDataSet* Foam::vtkPV3Foam::GetDataSetFromBlock
(
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);
    if (block)
    {
        return vtkDataSet::SafeDownCast(block->GetBlock(datasetNo));
    }

    return 0;
}


Foam::label Foam::vtkPV3Foam::GetNumberOfDataSets
(
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);
    if (block)
    {
        return block->GetNumberOfBlocks();
    }

    return 0;
}


void Foam::vtkPV3Foam::resetCounters()
{
    // Reset region ids and sizes
    selectInfoVolume_.reset();
    selectInfoPatches_.reset();
    selectInfoLagrangian_.reset();
    selectInfoCellZones_.reset();
    selectInfoFaceZones_.reset();
    selectInfoPointZones_.reset();
    selectInfoCellSets_.reset();
    selectInfoFaceSets_.reset();
    selectInfoPointSets_.reset();
}


bool Foam::vtkPV3Foam::setTime(const double& requestedTime)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::setTime(" << requestedTime << ")"
            << endl;
    }

    Time& runTime = dbPtr_();

    // Get times list
    instantList Times = runTime.times();

    bool found = false;
    int nearestIndex = Time::findClosestTimeIndex(Times, requestedTime);

    if (nearestIndex == -1)
    {
        nearestIndex = 0;
        found = false;
    }
    else
    {
        found = true;
    }

    // see what has changed
    if (timeIndex_ != nearestIndex)
    {
        timeIndex_ = nearestIndex;
        runTime.setTime(Times[nearestIndex], nearestIndex);

        // the fields change each time
        fieldsChanged_ = true;

        if (meshPtr_)
        {
            if (meshPtr_->readUpdate() != polyMesh::UNCHANGED)
            {
                meshChanged_ = true;
                // patches, zones etc might have changed
                UpdateInformation();
            }
        }
        else
        {
            meshChanged_ = true;
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::setTime() - selected time "
            << Times[nearestIndex].name() << " index=" << nearestIndex
            << " meshChanged=" << meshChanged_
            << " fieldsChanged=" << fieldsChanged_ << endl;
    }

    return found;
}


void Foam::vtkPV3Foam::updateSelectedRegions()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateSelectedRegions" << endl;
    }

    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    const label nSelect = arraySelection->GetNumberOfArrays();

    if (selectedRegions_.size() != nSelect)
    {
        selectedRegions_.setSize(nSelect);
        selectedRegions_ = 0;
        meshChanged_ = true;
    }

    selectedRegionDatasetIds_.setSize(nSelect);

    // Read the selected cell regions, zones, patches and add to region list
    forAll (selectedRegions_, regionId)
    {
        int setting = arraySelection->GetArraySetting(regionId);

        if (selectedRegions_[regionId] != setting)
        {
            selectedRegions_[regionId] = setting;
            meshChanged_ = true;
        }

        selectedRegionDatasetIds_[regionId] = -1;

        if (debug)
        {
            Info<< "  region[" << regionId << "] = "
                << selectedRegions_[regionId]
                << " : " << arraySelection->GetArrayName(regionId) << endl;
        }
    }
    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateSelectedRegions" << endl;
    }
}


Foam::stringList Foam::vtkPV3Foam::getSelectedArrayEntries
(
    vtkDataArraySelection* arraySelection,
    const bool firstWord
)
{
    stringList selections(arraySelection->GetNumberOfArrays());
    label nElem = 0;

    if (debug)
    {
        Info<< "available(";
        forAll (selections, elemI)
        {
            Info<< " \"" << arraySelection->GetArrayName(elemI) << "\"";
        }
        Info<< " )\n"
            << "selected(";
    }

    forAll (selections, elemI)
    {
        if (arraySelection->GetArraySetting(elemI))
        {
            if (firstWord)
            {
                selections[nElem] = getFirstWord
                (
                    arraySelection->GetArrayName(elemI)
                );
            }
            else
            {
                selections[nElem] = arraySelection->GetArrayName(elemI);
            }

            if (debug)
            {
                Info<< " " << selections[nElem];
            }

            ++nElem;
        }
    }

    if (debug)
    {
        Info<< " )" << endl;
    }

    selections.setSize(nElem);
    return selections;
}


Foam::stringList Foam::vtkPV3Foam::getSelectedArrayEntries
(
    vtkDataArraySelection* arraySelection,
    const selectionInfo& selector,
    const bool firstWord
)
{
    stringList selections(selector.size());
    label nElem = 0;

    if (debug)
    {
        Info<< "available(";
        for
        (
            int elemI = selector.start();
            elemI < selector.end();
            ++elemI
        )
        {
            Info<< " \"" << arraySelection->GetArrayName(elemI) << "\"";
        }

        Info<< " )\n"
            << "selected(";
    }

    for
    (
        int elemI = selector.start();
        elemI < selector.end();
        ++elemI
    )
    {
        if (arraySelection->GetArraySetting(elemI))
        {
            if (firstWord)
            {
                selections[nElem] = getFirstWord
                (
                    arraySelection->GetArrayName(elemI)
                );
            }
            else
            {
                selections[nElem] = arraySelection->GetArrayName(elemI);
            }

            if (debug)
            {
                Info<< " " << selections[nElem];
            }

            ++nElem;
        }
    }

    if (debug)
    {
        Info<< " )" << endl;
    }

    selections.setSize(nElem);
    return selections;
}


void Foam::vtkPV3Foam::setSelectedArrayEntries
(
    vtkDataArraySelection* arraySelection,
    const stringList& selections
)
{
    if (debug > 1)
    {
        Info<< "<beg> Foam::vtkPV3Foam::setSelectedArrayEntries" << endl;
    }
    const label nEntries = arraySelection->GetNumberOfArrays();

    // Reset all current entries to 'not selected'
    arraySelection->DisableAllArrays();

    // Loop through entries, setting values from selectedEntries
    forAll (selections, elemI)
    {
        if (debug > 1)
        {
            Info<< "selections[" << elemI << "] = " << selections[elemI]
                << endl;
        }

        for (label i=0; i<nEntries; i++)
        {
            string arrayName = arraySelection->GetArrayName(i);

            if (arrayName == selections[elemI])
            {
                if (debug > 1)
                {
                    Info<< "enabling array: " << arrayName << " Index = "
                        << i
                        << endl;
                }

                arraySelection->EnableArray(arrayName.c_str());
                break;
            }
        }
    }
    if (debug > 1)
    {
        Info<< "<end> Foam::vtkPV3Foam::setSelectedArrayEntries" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPV3Foam::vtkPV3Foam
(
    const char* const FileName,
    vtkPV3FoamReader* reader
)
:
    reader_(reader),
    dbPtr_(NULL),
    meshPtr_(NULL),
    selectInfoVolume_(VOLUME, "unzoned"),
    selectInfoPatches_(PATCHES, "patches"),
    selectInfoLagrangian_(LAGRANGIAN, "lagrangian"),
    selectInfoCellZones_(CELLZONE, "cellZone"),
    selectInfoFaceZones_(FACEZONE, "faceZone"),
    selectInfoPointZones_(POINTZONE, "pointZone"),
    selectInfoCellSets_(CELLSET, "cellSet"),
    selectInfoFaceSets_(FACESET, "faceSet"),
    selectInfoPointSets_(POINTSET, "pointSet"),
    patchTextActorsPtrs_(0),
    nMesh_(0),
    timeIndex_(-1),
    meshChanged_(true),
    fieldsChanged_(true)
{
    if (debug)
    {
        Info<< "Foam::vtkPV3Foam::vtkPV3Foam - " << FileName << endl;
        printMemory();
    }

    // avoid argList and get rootPath/caseName directly from the file
    fileName fullCasePath(fileName(FileName).path());

    if (!dir(fullCasePath))
    {
        return;
    }
    if (fullCasePath == ".")
    {
        fullCasePath = cwd();
    }

    // Set the case as an environment variable - some BCs might use this
    if (fullCasePath.name().find("processor", 0) == 0)
    {
        setEnv("FOAM_CASE", fullCasePath.path(), true);
    }
    else
    {
        setEnv("FOAM_CASE", fullCasePath, true);
    }

    if (debug)
    {
        Info<< "fullCasePath=" << fullCasePath << nl
            << "FOAM_CASE=" << getEnv("FOAM_CASE") << endl;
    }

    // Create time object
    dbPtr_.reset
    (
        new Time
        (
            Time::controlDictName,
            fileName(fullCasePath.path()),
            fileName(fullCasePath.name())
        )
    );

    dbPtr_().functionObjects().off();

    // Set initial cloud name
    // TODO - TEMPORARY MEASURE UNTIL CAN PROCESS MULTIPLE CLOUDS
    cloudName_ = "";

    UpdateInformation();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkPV3Foam::~vtkPV3Foam()
{
    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::~vtkPV3Foam" << endl;
    }

    if (meshPtr_)
    {
        delete meshPtr_;
        meshPtr_ = NULL;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::UpdateInformation()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::UpdateInformation"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "] TimeStep="
            << reader_->GetTimeStep() << endl;
    }

    resetCounters();

    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    stringList selectedEntries;
    // enable 'internalMesh' on the first call
    if (arraySelection->GetNumberOfArrays() == 0 && !meshPtr_)
    {
        selectedEntries.setSize(1);
        selectedEntries[0] = "internalMesh";
    }
    else
    {
        // preserve the currently selected values
        selectedEntries = getSelectedArrayEntries
        (
            arraySelection
        );
    }

    // Clear current region list/array
    arraySelection->RemoveAllArrays();

    // Update region array
    updateInformationInternalMesh();
    updateInformationLagrangian();
    updateInformationPatches();

    if (reader_->GetIncludeSets())
    {
        updateInformationSets();
    }

    if (reader_->GetIncludeZones())
    {
        updateInformationZones();
    }

    // restore the currently enabled values
    setSelectedArrayEntries
    (
        arraySelection,
        selectedEntries
    );

    if (meshChanged_)
    {
        fieldsChanged_ = true;
    }

    // Update volField array
    updateInformationFields<fvPatchField, volMesh>
    (
        reader_->GetVolFieldSelection()
    );

    // Update pointField array
    updateInformationFields<pointPatchField, pointMesh>
    (
        reader_->GetPointFieldSelection()
    );

    // Update lagrangian field array
    updateInformationLagrangianFields();

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::UpdateInformation" << endl;
    }

}


void Foam::vtkPV3Foam::Update
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        cout<< "<beg> Foam::vtkPV3Foam::Update" << nl
            <<"Update\n";
        output->Print(cout);

        cout<<"Internally:\n";
        output->Print(cout);

        cout<< " has " << output->GetNumberOfBlocks() << " blocks\n";
        printMemory();
    }

    // Set up region selection(s)
    updateSelectedRegions();

    // Update the Foam mesh
    updateFoamMesh();
    reader_->UpdateProgress(0.2);

    // Convert meshes
    convertMeshVolume(output);
    convertMeshLagrangian(output);
    convertMeshPatches(output);
    reader_->UpdateProgress(0.4);

    if (reader_->GetIncludeZones())
    {
        convertMeshCellZones(output);
        convertMeshFaceZones(output);
        convertMeshPointZones(output);
    }

    if (reader_->GetIncludeSets())
    {
        convertMeshCellSets(output);
        convertMeshFaceSets(output);
        convertMeshPointSets(output);
    }
    reader_->UpdateProgress(0.8);

    // Update fields
    updateVolFields(output);
    updatePointFields(output);
    updateLagrangianFields(output);
    reader_->UpdateProgress(1.0);

    if (debug)
    {
        Info<< "Number of data sets after update" << nl
            << "  VOLUME = "
            << GetNumberOfDataSets(output, selectInfoVolume_) << nl
            << "  PATCHES = "
            << GetNumberOfDataSets(output, selectInfoPatches_) << nl
            << "  LAGRANGIAN = "
            << GetNumberOfDataSets(output, selectInfoLagrangian_) << nl
            << "  CELLZONE = "
            << GetNumberOfDataSets(output, selectInfoCellZones_) << nl
            << "  FACEZONE = "
            << GetNumberOfDataSets(output, selectInfoFaceZones_) << nl
            << "  POINTZONE = "
            << GetNumberOfDataSets(output, selectInfoPointZones_) << nl
            << "  CELLSET = "
            << GetNumberOfDataSets(output, selectInfoCellSets_) << nl
            << "  FACESET = "
            << GetNumberOfDataSets(output, selectInfoFaceSets_) << nl
            << "  POINTSET = "
            << GetNumberOfDataSets(output, selectInfoPointSets_) << nl;

        // traverse blocks:
        cout<< "nBlocks = " << output->GetNumberOfBlocks() << "\n";
        cout<< "done Update\n";
        output->Print(cout);
        cout<< " has " << output->GetNumberOfBlocks() << " blocks\n";
        output->GetInformation()->Print(cout);

        cout<<"ShouldIReleaseData :" << output->ShouldIReleaseData() << "\n";
        printMemory();
    }

    meshChanged_ = fieldsChanged_ = false;
}


double* Foam::vtkPV3Foam::findTimes(int& nTimeSteps)
{
    int nTimes = 0;
    double* tsteps = NULL;

    if (dbPtr_.valid())
    {
        Time& runTime = dbPtr_();
        instantList timeLst = runTime.times();

        // always skip "constant" time, unless there are no other times
        nTimes = timeLst.size();
        label timeI = 0;

        if (nTimes > 1)
        {
            timeI = 1;
            --nTimes;
        }

        if (nTimes)
        {
            tsteps = new double[nTimes];
            for (label stepI = 0; stepI < nTimes; ++stepI, ++timeI)
            {
                tsteps[stepI] = timeLst[timeI].value();
            }
        }
    }
    else
    {
        if (debug)
        {
            cout<< "no valid dbPtr:\n";
        }
    }

    // vector length returned via the parameter
    nTimeSteps = nTimes;

    return tsteps;
}


void Foam::vtkPV3Foam::addPatchNames(vtkRenderer* renderer)
{
    // Remove any patch names previously added to the renderer
    removePatchNames(renderer);

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::addPatchNames" << endl;
    }

    const polyBoundaryMesh& pbMesh = meshPtr_->boundaryMesh();

    const selectionInfo& selector = selectInfoPatches_;

    // the currently selected patches, strip off any suffix
    const stringList selectedPatches = getSelectedArrayEntries
    (
        reader_->GetRegionSelection(),
        selector,
        true
    );

    if (debug)
    {
        Info<<"... add patches: " << selectedPatches << endl;
    }

    // Find the total number of zones
    // Each zone will take the patch name

    // Number of zones per patch ... zero zones should be skipped
    labelList nZones(pbMesh.size(), 0);

    // Per global zone number the average face centre position
    DynamicList<point> zoneCentre(pbMesh.size());

    if (debug)
    {
        Info<< "... determining patch zones" << endl;
    }

    // Loop through all patches to determine zones, and centre of each zone
    forAll(pbMesh, patchI)
    {
        const polyPatch& pp = pbMesh[patchI];

        // Only include the patch if it is selected
        bool isSelected = false;
        forAll(selectedPatches, elemI)
        {
            if (pp.name() == selectedPatches[elemI])
            {
                isSelected = true;
                break;
            }
        }

        if (isSelected)
        {
            const labelListList& edgeFaces = pp.edgeFaces();
            const vectorField& n = pp.faceNormals();

            boolList featEdge(pp.nEdges(), false);

            forAll(edgeFaces, edgeI)
            {
                const labelList& eFaces = edgeFaces[edgeI];

                if (eFaces.size() != 2)
                {
                    featEdge[edgeI] = true;
                }
                else if (mag(n[eFaces[0]] & n[eFaces[1]]) < 0.5)
                {
                    featEdge[edgeI] = true;
                }
            }

            // Do topological analysis of patch. Determine disconnected regions
            patchZones pZones(pp, featEdge);

            nZones[patchI] = pZones.nZones();

            labelList zoneNFaces(pZones.nZones(), 0);

            // Save start of information for current patch
            label patchStart = zoneCentre.size();

            // Create storage for additional zone centres
            forAll(zoneNFaces, zoneI)
            {
                zoneCentre.append(vector::zero);
            }

            // Do averaging per individual zone

            forAll(pp, faceI)
            {
                label zoneI = pZones[faceI];
                zoneCentre[patchStart+zoneI] += pp[faceI].centre(pp.points());
                zoneNFaces[zoneI]++;
            }

            for (label i=0; i<nZones[patchI]; i++)
            {
                zoneCentre[patchStart + i] /= zoneNFaces[i];
            }
        }
    }
    zoneCentre.shrink();

    if (debug)
    {
        Info<< "patch zone centres = " << zoneCentre << nl
            << "zones per patch = " << nZones << endl;
    }

    // Set the size of the patch labels to max number of zones
    patchTextActorsPtrs_.setSize(zoneCentre.size());

    if (debug)
    {
        Info<< "constructing patch labels" << endl;
    }

    label globalZoneI = 0;
    forAll(pbMesh, patchI)
    {
        const polyPatch& pp = pbMesh[patchI];

        // Only selected patches will have a non-zero number of zones
        for (label i=0; i<nZones[patchI]; i++)
        {
            if (debug)
            {
                Info<< "patch name = " << pp.name() << nl
                    << "anchor = " << zoneCentre[globalZoneI] << nl
                    << "globalZoneI = " << globalZoneI << endl;
            }

            vtkTextActor* txt = vtkTextActor::New();

            txt->SetInput(pp.name().c_str());

            // Set text properties
            vtkTextProperty* tprop = txt->GetTextProperty();
            tprop->SetFontFamilyToArial();
            tprop->BoldOff();
            tprop->ShadowOff();
            tprop->SetLineSpacing(1.0);
            tprop->SetFontSize(12);
            tprop->SetColor(1.0, 0.0, 0.0);
            tprop->SetJustificationToCentered();

            // Set text to use 3-D world co-ordinates
            txt->GetPositionCoordinate()->SetCoordinateSystemToWorld();

            txt->GetPositionCoordinate()->SetValue
            (
                zoneCentre[globalZoneI].x(),
                zoneCentre[globalZoneI].y(),
                zoneCentre[globalZoneI].z()
            );

            // Add text to each renderer
            renderer->AddViewProp(txt);

            // Maintain a list of text labels added so that they can be
            // removed later
            patchTextActorsPtrs_[globalZoneI] = txt;

            globalZoneI++;
        }
    }

    // Resize the patch names list to the actual number of patch names added
    patchTextActorsPtrs_.setSize(globalZoneI);

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::addPatchNames)" << endl;
    }
}


void Foam::vtkPV3Foam::removePatchNames(vtkRenderer* renderer)
{
    if (debug)
    {
        Info<< "removePatchNames()" << endl;
    }

    forAll(patchTextActorsPtrs_, patchI)
    {
        renderer->RemoveViewProp(patchTextActorsPtrs_[patchI]);
        patchTextActorsPtrs_[patchI]->Delete();
    }
    patchTextActorsPtrs_.setSize(0);
}


void Foam::vtkPV3Foam::PrintSelf(ostream& os, vtkIndent indent) const
{
    os  << indent << "Number of meshes: " << nMesh_ << "\n";
    os  << indent << "Number of nodes: "
        << (meshPtr_ ? meshPtr_->nPoints() : 0) << "\n";

    os  << indent << "Number of cells: "
        << (meshPtr_ ? meshPtr_->nCells() : 0) << "\n";

    os  << indent << "Number of available time steps: "
        << (dbPtr_.valid() ? dbPtr_().times().size() : 0) << endl;
}


// parse these bits of info from /proc/meminfo (Linux)
//
// MemTotal:      2062660 kB
// MemFree:       1124400 kB
//
// used = MemTotal - MemFree is what the free(1) uses.
//
void Foam::vtkPV3Foam::printMemory()
{
    const char* meminfo = "/proc/meminfo";

    if (exists(meminfo))
    {
        IFstream is(meminfo);
        label memTotal = 0;
        label memFree = 0;

        string line;

        while (is.getLine(line).good())
        {
            char tag[32];
            int value;

            if (sscanf(line.c_str(), "%30s %d", tag, &value) == 2)
            {
                if (!strcmp(tag, "MemTotal:"))
                {
                    memTotal = value;
                }
                else if (!strcmp(tag, "MemFree:"))
                {
                    memFree = value;
                }
            }
        }

        Info << "memUsed: " << (memTotal - memFree) << " kB" << endl;
    }
    else
    {
        Info << "Has no /proc/meminfo - Not a Linux-machine" << endl;
    }
}

// ************************************************************************* //
