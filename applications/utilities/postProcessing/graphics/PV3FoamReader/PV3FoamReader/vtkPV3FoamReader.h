/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPV3FoamReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPV3FoamReader - reads a dataset in OpenFOAM format
// .SECTION Description
// vtkPV3FoamReader creates an multiblock dataset.
// It uses the OpenFOAM infrastructure (fvMesh, etc) to
// handle mesh and field data.

#ifndef __vtkPV3FoamReader_h
#define __vtkPV3FoamReader_h

#include "vtkMultiBlockDataSetAlgorithm.h"

// Foam forward declarations
namespace Foam
{
    class vtkPV3Foam;
}

// VTK forward declarations
class vtkUnstructuredGrid;
class vtkPoints;
class vtkIntArray;
class vtkFloatArray;
class vtkDoubleArray;
class vtkDataArraySelection;
class vtkCallbackCommand;


class VTK_IO_EXPORT vtkPV3FoamReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:

    static vtkPV3FoamReader* New();

    vtkTypeRevisionMacro
    (
        vtkPV3FoamReader,
        vtkMultiBlockDataSetAlgorithm
    );

    void PrintSelf
    (
        ostream& os,
        vtkIndent indent
    );

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Time control
    // Set/Get the timestep and the timestep range
    vtkSetMacro(TimeStep, int);
    vtkGetMacro(TimeStep, int);
    vtkSetVector2Macro(TimeStepRange, int);
    vtkGetVector2Macro(TimeStepRange, int);

    // GUI update control
    vtkSetMacro(UpdateGUI, int);
    vtkGetMacro(UpdateGUI, int);

    // FOAM mesh caching control
    vtkSetMacro(CacheMesh, int);
    vtkGetMacro(CacheMesh, int);

    // FOAM extrapolate internal values onto the walls
    vtkSetMacro(ExtrapolateWalls, int);
    vtkGetMacro(ExtrapolateWalls, int);

    // FOAM read sets control
    vtkSetMacro(IncludeSets, int);
    vtkGetMacro(IncludeSets, int);

    // FOAM read zones control
    vtkSetMacro(IncludeZones, int);
    vtkGetMacro(IncludeZones, int);

    // FOAM display patch names control
    vtkSetMacro(ShowPatchNames, int);
    vtkGetMacro(ShowPatchNames, int);

    // Region selection list control
    vtkDataArraySelection* GetRegionSelection();
    int GetNumberOfRegionArrays();
    const char* GetRegionArrayName(int index);
    int GetRegionArrayStatus(const char* name);
    void SetRegionArrayStatus(const char* name, int status);

    // volField selection list control
    vtkDataArraySelection* GetVolFieldSelection();
    int GetNumberOfVolFieldArrays();
    const char* GetVolFieldArrayName(int index);
    int GetVolFieldArrayStatus(const char* name);
    void SetVolFieldArrayStatus(const char* name, int status);

    // pointField selection list control
    vtkDataArraySelection* GetPointFieldSelection();
    int GetNumberOfPointFieldArrays();
    int GetPointFieldArrayStatus(const char* name);
    void SetPointFieldArrayStatus(const char* name, int status);
    const char* GetPointFieldArrayName(int index);

    // lagrangianField selection list control
    vtkDataArraySelection* GetLagrangianFieldSelection();
    int GetNumberOfLagrangianFieldArrays();
    int GetLagrangianFieldArrayStatus(const char* name);
    void SetLagrangianFieldArrayStatus(const char* name, int status);
    const char* GetLagrangianFieldArrayName(int index);

    // Callback registered with the SelectionObserver
    // for all the selection lists
    static void SelectionModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );

    void SelectionModified();


protected:

    vtkPV3FoamReader();
    ~vtkPV3FoamReader();

    char* FileName;

    virtual int RequestData
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    virtual int RequestInformation
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    // The observer to modify this object when the array selections
    // are modified
    vtkCallbackCommand* SelectionObserver;


private:

    vtkPV3FoamReader(const vtkPV3FoamReader&);  // Not implemented.
    void operator=(const vtkPV3FoamReader&);  // Not implemented.

    //- Add patch names to the view
    void addPatchNamesToView();

    //- Remove patch names from the view
    void removePatchNamesFromView();

    int TimeStep;
    int TimeStepRange[2];

    int CacheMesh;

    int ExtrapolateWalls;
    int IncludeSets;
    int IncludeZones;
    int ShowPatchNames;

    int UpdateGUI;
    int UpdateGUIOld;

    vtkDataArraySelection* RegionSelection;
    vtkDataArraySelection* VolFieldSelection;
    vtkDataArraySelection* PointFieldSelection;
    vtkDataArraySelection* LagrangianFieldSelection;

    //- Access to the output port1
    vtkMultiBlockDataSet* output1_;

    //BTX
    Foam::vtkPV3Foam* foamData_;
    //ETX
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
