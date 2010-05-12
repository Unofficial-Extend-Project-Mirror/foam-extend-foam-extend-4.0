/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPV3FoamReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPV3FoamReader.h"

#include "pqApplicationCore.h"
#include "pqRenderView.h"
#include "pqServerManagerModel.h"

// VTK includes
#include "vtkCallbackCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkSMRenderViewProxy.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"

// Foam includes
#include "vtkPV3Foam.H"

vtkCxxRevisionMacro(vtkPV3FoamReader, "$Revision: 1.5$");
vtkStandardNewMacro(vtkPV3FoamReader);


vtkPV3FoamReader::vtkPV3FoamReader()
{
    Debug = 0;
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName  = NULL;
    foamData_ = NULL;

    output1_ = NULL;

    TimeStep = 0;
    TimeStepRange[0] = 0;
    TimeStepRange[1] = 0;

    CacheMesh = 0;

    ExtrapolateWalls = 0;
    IncludeSets = 0;
    IncludeZones = 0;
    ShowPatchNames = 0;

    UpdateGUI = 1;
    UpdateGUIOld = 1;

    RegionSelection = vtkDataArraySelection::New();
    VolFieldSelection = vtkDataArraySelection::New();
    PointFieldSelection = vtkDataArraySelection::New();
    LagrangianFieldSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPV3FoamReader::SelectionModifiedCallback
    );
    SelectionObserver->SetClientData(this);

    RegionSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    VolFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    PointFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    LagrangianFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
}


vtkPV3FoamReader::~vtkPV3FoamReader()
{
    vtkDebugMacro(<<"Deconstructor");

    if (foamData_)
    {
        delete foamData_;
    }

    if (FileName)
    {
        delete [] FileName;
    }

    RegionSelection->RemoveObserver(this->SelectionObserver);
    VolFieldSelection->RemoveObserver(this->SelectionObserver);
    PointFieldSelection->RemoveObserver(this->SelectionObserver);
    LagrangianFieldSelection->RemoveObserver(this->SelectionObserver);

    SelectionObserver->Delete();

    RegionSelection->Delete();
    VolFieldSelection->Delete();
    PointFieldSelection->Delete();
    LagrangianFieldSelection->Delete();
}


// Do everything except set the output info
int vtkPV3FoamReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestInformation");

    if (Foam::vtkPV3Foam::debug)
    {
        cout<<"REQUEST_INFORMATION\n";
    }

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    if (!foamData_)
    {
        vtkDebugMacro("RequestInformation: creating foamData_");
        vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
        (
            outInfo->Get(vtkMultiBlockDataSet::DATA_OBJECT())
        );

        if (Foam::vtkPV3Foam::debug)
        {
            cout<< "constructed vtkPV3Foam with output: ";
            output->Print(cout);
        }

        foamData_ = new Foam::vtkPV3Foam(FileName, this);
    }
    else
    {
        foamData_->UpdateInformation();
    }

    int nTimeSteps = 0;
    double* timeSteps = foamData_->findTimes(nTimeSteps);

    outInfo->Set
    (
        vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
        timeSteps,
        nTimeSteps
    );

    double timeRange[2];
    if (nTimeSteps)
    {
        timeRange[0] = timeSteps[0];
        timeRange[1] = timeSteps[nTimeSteps-1];

        if (Foam::vtkPV3Foam::debug > 1)
        {
            cout<<"nTimeSteps " << nTimeSteps << "\n";
            cout<<"timeRange " << timeRange[0] << " to " << timeRange[1] << "\n";

            for (int i = 0; i < nTimeSteps; ++i)
            {
                cout<< "step[" << i << "] = " << timeSteps[i] << "\n";
            }
        }

        outInfo->Set
        (
            vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
            timeRange,
            2
        );
    }

    delete timeSteps;

    return 1;
}


// Set the output info
int vtkPV3FoamReader::RequestData
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestData");

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    if (Foam::vtkPV3Foam::debug)
    {
        int nInfo = outputVector->GetNumberOfInformationObjects();
        cout<<"requestData with " << nInfo << " items\n";

        for (int i = 0; i < nInfo; ++i)
        {
            vtkInformation *info = outputVector->GetInformationObject(i);
            info->Print(cout);
        }
    }

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
    (
        outInfo->Get(vtkMultiBlockDataSet::DATA_OBJECT())
    );

    if (Foam::vtkPV3Foam::debug)
    {
        vtkInformation* outputInfo = this->GetOutputPortInformation(0);
        outputInfo->Print(cout);

        vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
        (
            outputInfo->Get(vtkMultiBlockDataSet::DATA_OBJECT())
        );
        if (output)
        {
            output->Print(cout);
        }
        else
        {
            cout<< "no output\n";
        }

        vtkInformation* execInfo = this->GetExecutive()->GetOutputInformation(0);
        execInfo->Print(cout);

        outInfo->Print(cout);

        vtkMultiBlockDataSet* dobj = vtkMultiBlockDataSet::SafeDownCast
        (
            outInfo->Get(vtkMultiBlockDataSet::DATA_OBJECT())
        );
        if (dobj)
        {
            dobj->Print(cout);

            vtkInformation* dobjInfo = dobj->GetInformation();
            dobjInfo->Print(cout);
        }
        else
        {
            cout<< "no data_object\n";
        }
    }

    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
        if (Foam::vtkPV3Foam::debug)
        {
            cout<<"Has UPDATE_TIME_STEPS\n";
            cout<<"output->GetNumberOfBlocks() = "
                << output->GetNumberOfBlocks() << "\n";
        }

        // Get the requested time step.
        // We only supprt requests of a single time step
        int nRequestedTimeSteps = outInfo->Length
        (
            vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()
        );
        if (nRequestedTimeSteps >= 1)
        {
            double *requestedTimeSteps = outInfo->Get
            (
                vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()
            );

            foamData_->setTime(requestedTimeSteps[0]);
        }
    }

    if
    (
        (UpdateGUIOld == GetUpdateGUI())
     || (output->GetNumberOfBlocks() == 0)
    )
    {
        foamData_->Update(output);

        if (ShowPatchNames == 1)
        {
            addPatchNamesToView();
        }
        else
        {
            removePatchNamesFromView();
        }
    }
    UpdateGUIOld = GetUpdateGUI();

    return 1;
}


void vtkPV3FoamReader::addPatchNamesToView()
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    if (appCore==NULL)
    {
        return;
    }

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();

    // Get all the pqRenderView instances
    QList<pqRenderView*> renderViews = smModel->findItems<pqRenderView*>();

    for (int viewI=0; viewI<renderViews.size(); viewI++)
    {
        foamData_->addPatchNames
        (
            renderViews[viewI]->getRenderViewProxy()->GetRenderer()
        );
    }
}


void vtkPV3FoamReader::removePatchNamesFromView()
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    if (appCore==NULL)
    {
        return;
    }

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();

    // Get all the pqRenderView instances
    QList<pqRenderView*> renderViews = smModel->findItems<pqRenderView*>();

    for (int viewI=0; viewI<renderViews.size(); viewI++)
    {
        foamData_->removePatchNames
        (
            renderViews[viewI]->getRenderViewProxy()->GetRenderer()
        );
    }
}


void vtkPV3FoamReader::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os<< indent << "File name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";

    foamData_->PrintSelf(os, indent);

    os<< indent << "Time step range: "
      << this->TimeStepRange[0] << " - " << this->TimeStepRange[1]
      << "\n";
    os<< indent << "Time step: " << this->TimeStep << endl;
}


// ----------------------------------------------------------------------
// Region selection list control

vtkDataArraySelection* vtkPV3FoamReader::GetRegionSelection()
{
    vtkDebugMacro(<<"GetRegionSelection");
    return RegionSelection;
}


int vtkPV3FoamReader::GetNumberOfRegionArrays()
{
    vtkDebugMacro(<<"GetNumberOfRegionArrays");
    return RegionSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetRegionArrayName(int index)
{
    vtkDebugMacro(<<"GetRegionArrayName");
    return RegionSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetRegionArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetRegionArrayStatus");
    return RegionSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetRegionArrayStatus(const char* name, int status)
{
    vtkDebugMacro(<<"SetRegionArrayStatus");

    if (status)
    {
        RegionSelection->EnableArray(name);
    }
    else
    {
        RegionSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// volField selection list control

vtkDataArraySelection* vtkPV3FoamReader::GetVolFieldSelection()
{
    vtkDebugMacro(<<"GetVolFieldSelection");
    return VolFieldSelection;
}


int vtkPV3FoamReader::GetNumberOfVolFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfVolFieldArrays");
    return VolFieldSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetVolFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetVolFieldArrayName");
    return VolFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetVolFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetVolFieldArrayStatus");
    return VolFieldSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetVolFieldArrayStatus(const char* name, int status)
{
    vtkDebugMacro(<<"SetVolFieldArrayStatus");
    if (status)
    {
        VolFieldSelection->EnableArray(name);
    }
    else
    {
        VolFieldSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// pointField selection list control

vtkDataArraySelection* vtkPV3FoamReader::GetPointFieldSelection()
{
    vtkDebugMacro(<<"GetPointFieldSelection");
    return PointFieldSelection;
}


int vtkPV3FoamReader::GetNumberOfPointFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfPointFieldArrays");
    return PointFieldSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetPointFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetPointFieldArrayName");
    return PointFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetPointFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetPointFieldArrayStatus");
    return PointFieldSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetPointFieldArrayStatus(const char* name, int status)
{
    vtkDebugMacro(<<"SetPointFieldArrayStatus");
    if (status)
    {
        PointFieldSelection->EnableArray(name);
    }
    else
    {
        PointFieldSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// lagrangianField selection list control

vtkDataArraySelection* vtkPV3FoamReader::GetLagrangianFieldSelection()
{
    vtkDebugMacro(<<"GetLagrangianFieldSelection");
    return LagrangianFieldSelection;
}


int vtkPV3FoamReader::GetNumberOfLagrangianFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfLagrangianFieldArrays");
    return LagrangianFieldSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetLagrangianFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayName");
    return LagrangianFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetLagrangianFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayStatus");
    return LagrangianFieldSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetLagrangianFieldArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetLagrangianFieldArrayStatus");
    if (status)
    {
        LagrangianFieldSelection->EnableArray(name);
    }
    else
    {
        LagrangianFieldSelection->DisableArray(name);
    }
}

// ----------------------------------------------------------------------

void vtkPV3FoamReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkPV3FoamReader*>(clientdata)->SelectionModified();
}


void vtkPV3FoamReader::SelectionModified()
{
    vtkDebugMacro(<<"SelectionModified");
    Modified();
}

// ************************************************************************* //
