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
    Foam::vtkPV3Foam

Description
    Provides a reader interface for OpenFOAM to VTK interaction.

SourceFiles
    vtkPV3Foam.C
    vtkPV3Foam.H
    vtkPV3FoamI.H
    vtkPV3FoamAddLagrangianMesh.C
    vtkPV3FoamAddPatchMesh.C
    vtkPV3FoamAddSetMesh.C
    vtkPV3FoamAddToSelection.H
    vtkPV3FoamAddVolumeMesh.C
    vtkPV3FoamAddZoneMesh.C
    vtkPV3FoamConvertFaceField.H
    vtkPV3FoamConvertLagrangianFields.H
    vtkPV3FoamConvertMesh.C
    vtkPV3FoamConvertPatchFaceField.H
    vtkPV3FoamConvertPatchPointField.H
    vtkPV3FoamConvertPointFields.H
    vtkPV3FoamConvertVolFields.H
    vtkPV3FoamInsertNextPoint.H
    vtkPV3FoamUpdate.C
    vtkPV3FoamUpdateInformation.C
    vtkPV3FoamUpdateInformationFields.H

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3Foam_H
#define vtkPV3Foam_H

#include "className.H"
#include "fileName.H"
#include "volPointInterpolation.H"
#include "stringList.H"
#include "wordList.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkPoints;
class vtkPV3FoamReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkUnstructuredGrid;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Foam class forward declarations
class argList;
class Time;
class fvMesh;
class IOobjectList;
class polyPatch;
class faceSet;
class pointSet;

template<class Type>
class IOField;

template<class Type>
class List;

/*---------------------------------------------------------------------------*\
                        Class vtkPV3Foam Declaration
\*---------------------------------------------------------------------------*/

class vtkPV3Foam
{
public:

    // Public data

        //- bookkeeping for the GUI checklists and the multi-block organization
        class selectionInfo
        {
            int block_;
            const char *name_;
            int start_;
            int size_;

        public:

            selectionInfo(const int blockNo, const char *name)
            :
                block_(blockNo),
                name_(name),
                start_(-1),
                size_(0)
            {}

            //- Return the block holding these datasets
            int block() const
            {
                return block_;
            }

            const char* name() const
            {
                return name_;
            }

            int start() const
            {
                return start_;
            }

            int end() const
            {
                return start_ + size_;
            }

            int size() const
            {
                return size_;
            }

            void reset()
            {
                start_ = -1;
                size_ = 0;
            }

            //- Assign new start and reset the size
            void operator=(const int i)
            {
                start_ = i;
                size_ = 0;
            }

            //- Increment the size
            void operator+=(const int n)
            {
                size_ += n;
            }
        };


private:

    // Private data

        //BTX
        //- Indices for datasets in vtkMultiBlockDataSet
        enum
        {
            VOLUME      = 0, // internal mesh
            PATCHES     = 1, // patches
            LAGRANGIAN  = 2,
            CELLZONE    = 3,
            FACEZONE    = 4,
            POINTZONE   = 5,
            CELLSET     = 6,
            FACESET     = 7,
            POINTSET    = 8
        };
        //ETX

        //- Access to the controlling vtkPV3FoamReader
        vtkPV3FoamReader* reader_;

        //- Foam time control
        autoPtr<Time> dbPtr_;

        //- Foam mesh
        fvMesh* meshPtr_;

        //- First instance and size of various regions
        selectionInfo selectInfoVolume_;
        selectionInfo selectInfoPatches_;
        selectionInfo selectInfoLagrangian_;
        selectionInfo selectInfoCellZones_;
        selectionInfo selectInfoFaceZones_;
        selectionInfo selectInfoPointZones_;
        selectionInfo selectInfoCellSets_;
        selectionInfo selectInfoFaceSets_;
        selectionInfo selectInfoPointSets_;

        //- Selected regions
        //  [0] = internal mesh, then lagrangian, patches, zones, sets
        boolList selectedRegions_;

        //- Selected regions indices in each respective block
        labelList selectedRegionDatasetIds_;

        //- Labels of cell-centres used as additional points when decomposing
        // polyhedra
        labelList addPointCellLabels_;

        //- Label of original cell for decomposed cells
        //  - internal mesh
        labelList superCells_;

        //- Label of original cell for decomposed cells
        //  - cellZone meshes
        List<labelList> zoneSuperCells_;

        //- Label of original cell for decomposed cells
        //  - cellSet meshes
        List<labelList> csetSuperCells_;

        //- List of patch names
        List<vtkTextActor*> patchTextActorsPtrs_;


        // Dataset sizes

            //- Number of meshes
            // TODO - for info only - only set up to process ONE mesh
            int nMesh_;

            //- Cloud name to be processed
            // TODO - currently only set up to process ONE cloud
            word cloudName_;

            //- The time index
            int timeIndex_;

            //- Track changes in mesh geometry
            bool meshChanged_;

            //- Track changes in fields
            bool fieldsChanged_;


    // Private Member Functions

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static void AddToBlock
        (
            vtkMultiBlockDataSet* output,
            const selectionInfo&,
            const label datasetNo,
            vtkDataSet* dataset,
            const string& blockName = string::null
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static vtkDataSet* GetDataSetFromBlock
        (
            vtkMultiBlockDataSet* output,
            const selectionInfo&,
            const label datasetNo
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static label GetNumberOfDataSets
        (
            vtkMultiBlockDataSet* output,
            const selectionInfo&
        );

        //- Reset data counters
        void resetCounters();

        // Update information helper functions

            //- Update the regions selected in the GUI
            void updateSelectedRegions();

            //- Internal mesh info
            void updateInformationInternalMesh();

            //- Lagrangian info
            void updateInformationLagrangian();

            //- Patch info
            void updateInformationPatches();

            //- Set info
            void updateInformationSets();

            //- Read zone names for zoneType from file
            wordList readZoneNames(const word& zoneType);

            //- Zone info
            void updateInformationZones();

            //- Add to paraview array selection
            template<class Type>
            label addToSelection
            (
                vtkDataArraySelection *arraySelection,
                const IOobjectList&,
                const string& suffix = ""
            );

            //- Field info
            template<template<class> class patchType, class meshType>
            void updateInformationFields
            (
                vtkDataArraySelection *fieldSelection
            );

            //- Lagrangian field info
            void updateInformationLagrangianFields();


        // Update helper functions

            //- Foam mesh
            void updateFoamMesh();

            //- Volume fields
            void updateVolFields(vtkMultiBlockDataSet* output);

            //- Point fields
            void updatePointFields(vtkMultiBlockDataSet* output);

            //- Lagrangian fields
            void updateLagrangianFields(vtkMultiBlockDataSet* output);


        // Mesh conversion functions

            //- Volume mesh
            void convertMeshVolume(vtkMultiBlockDataSet* output);

            //- Lagrangian mesh
            void convertMeshLagrangian(vtkMultiBlockDataSet* output);

            //- Patch meshes
            void convertMeshPatches(vtkMultiBlockDataSet* output);

            //- Cell zone meshes
            void convertMeshCellZones(vtkMultiBlockDataSet* output);

            //- Face zone meshes
            void convertMeshFaceZones(vtkMultiBlockDataSet* output);

            //- Point zone meshes
            void convertMeshPointZones(vtkMultiBlockDataSet* output);

            //- Cell set meshes
            void convertMeshCellSets(vtkMultiBlockDataSet* output);

            //- Face set meshes
            void convertMeshFaceSets(vtkMultiBlockDataSet* output);

            //- Point set meshes
            void convertMeshPointSets(vtkMultiBlockDataSet* output);


        // Add mesh functions

            //- Add internal mesh/cell set meshes
            void addVolumeMesh
            (
                const fvMesh&,
                vtkUnstructuredGrid*,
                labelList& superCells
            );

            //- Add Lagrangian mesh
            void addLagrangianMesh
            (
                const fvMesh&,
                vtkPolyData*
            );

            //- Add patch mesh
            void addPatchMesh
            (
                const polyPatch&,
                vtkPolyData*
            );

            //- Add face zone mesh
            void addFaceZoneMesh
            (
                const fvMesh&,
                const labelList& faceLabels,
                vtkPolyData*
            );

            //- Add point zone
            void addPointZoneMesh
            (
                const fvMesh&,
                const labelList& pointLabels,
                vtkPolyData*
            );

            //- Add cell set mesh
            void addCellSetMesh
            (
                const fvMesh&,
                vtkUnstructuredGrid*
            );

            //- Add face set mesh
            void addFaceSetMesh
            (
                const fvMesh&,
                const faceSet&,
                vtkPolyData*
            );

            //- Add point mesh
            void addPointSetMesh
            (
                const fvMesh&,
                const pointSet&,
                vtkPolyData*
            );


        //- Add the fields in the selected time directory to the selection
        //  lists
        template<class GeoField>
        label addObjectsToSelection
        (
            vtkDataArraySelection* fieldSelection,
            const IOobjectList& objects,
            const string& suffix = string::null
        );


        // Convert Foam fields

            //- Volume fields - all types
            template<class Type>
            void convertVolFields
            (
                const fvMesh&,
                const volPointInterpolation& pInterp,
                const PtrList<PrimitivePatchInterpolation<primitivePatch> >&,
                const IOobjectList& objects,
                vtkDataArraySelection* fieldSelection,
                vtkMultiBlockDataSet* output
            );

            //- Volume field - all types except scalar
            template<class Type>
            void convertVolField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo,
                labelList& superCells
            );

            //- Patch field
            template<class Type>
            void convertPatchFaceField
            (
                const word& name,
                const Field<Type>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo
            );

            //- face set/zone field
            template<class Type>
            void convertFaceField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo,
                const fvMesh&,
                const labelList& faceLabels
            );

            //- face set/zone field
            template<class Type>
            void convertFaceField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo,
                const fvMesh&,
                const faceSet&
            );

            //- Lagrangian fields - all types
            template<class Type>
            void convertLagrangianFields
            (
                const fvMesh&,
                const IOobjectList& objects,
                vtkDataArraySelection *fieldSelection,
                vtkMultiBlockDataSet* output
            );

            //- Lagrangian field - all types except scalar
            template<class Type>
            void convertLagrangianField
            (
                const IOField<Type>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo
            );

            //- Point fields - all types
            template<class Type>
            void convertPointFields
            (
                const fvMesh& mesh,
                const IOobjectList& objects,
                vtkDataArraySelection *fieldSelection,
                vtkMultiBlockDataSet* output
            );

            //- Point fields - all types except scalar
            template<class Type>
            void convertPointField
            (
                const GeometricField<Type, pointPatchField, pointMesh>& ptf,
                const GeometricField<Type, fvPatchField, volMesh>& tf,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo
            );

            //- Patch point field - all types except scalar
            template<class Type>
            void convertPatchPointField
            (
                const word& name,
                const Field<Type>& tf,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo
            );


       // GUI selection helper functions

            //- Extract up to the first non-word characters
            inline static word getFirstWord(const char*);

            //- Store the current selection(s)
            static stringList getSelectedArrayEntries
            (
                vtkDataArraySelection* arraySelection,
                const bool firstWord = false
            );

            //- Store the current selection(s) for a sub-selection
            static stringList getSelectedArrayEntries
            (
                vtkDataArraySelection* arraySelection,
                const selectionInfo&,
                const bool firstWord = false
            );

            //- Set selection(s)
            static void setSelectedArrayEntries
            (
                vtkDataArraySelection* arraySelection,
                const stringList& selectedEntries
            );


        //- Disallow default bitwise copy construct
        vtkPV3Foam(const vtkPV3Foam&);

        //- Disallow default bitwise assignment
        void operator=(const vtkPV3Foam&);


public:

    //- Static data members

        ClassName("vtkPV3Foam");


    // Constructors

        //- Construct from components
        vtkPV3Foam
        (
            const char* const FileName,
            vtkPV3FoamReader* reader
        );


    //- Destructor

        ~vtkPV3Foam();


    // Member Functions

        //- Update
        void UpdateInformation();

        void Update(vtkMultiBlockDataSet* output);

        //- Allocate and return a list of selected times
        //  returns the count via the parameter
        double* findTimes(int& nTimeSteps);

        //- Add patch names to the display
        void addPatchNames(vtkRenderer* renderer);

        //- Remove patch names from the display
        void removePatchNames(vtkRenderer* renderer);

        //- set the runTime to the requested time
        //  sets to "constant" on error
        bool setTime(const double& requestedTime);


     // Access

        //- Debug information
        void PrintSelf(ostream&, vtkIndent) const;

        //- Simple memory used debugging information
        static void printMemory();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#    include "vtkPV3FoamI.H"

#endif

// ************************************************************************* //
