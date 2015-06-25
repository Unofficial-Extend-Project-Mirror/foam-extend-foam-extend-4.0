/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "workflowControls.H"
#include "polyMeshGen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

const std::map<word, label> workflowControls::workflowSteps_ =
    populateWorkflowSteps();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool workflowControls::restartRequested() const
{
    const dictionary& meshDict =
        mesh_.returnTime().lookupObject<dictionary>("meshDict");

    if
    (
        meshDict.found("workflowControls") &&
        meshDict.isDict("workflowControls")
    )
    {
        const dictionary& workflowControls =
            meshDict.subDict("workflowControls");

        if( workflowControls.found("restartFromLatestStep") )
        {
            const bool restart =
                readBool(workflowControls.lookup("restartFromLatestStep"));

            return restart;
        }
    }

    return false;
}

void workflowControls::setStepCompleted() const
{
    if( mesh_.metaData().found("lastStep") )
    {
        mesh_.metaData().set("lastStep", currentStep_);
    }
    else
    {
        mesh_.metaData().add("lastStep", currentStep_);
    }

    DynList<word> completedSteps;
    if( mesh_.metaData().found("completedSteps") )
        completedSteps = wordList(mesh_.metaData().lookup("completedSteps"));

    completedSteps.append(currentStep_);

    if( mesh_.metaData().found("completedSteps") )
    {
        mesh_.metaData().set("completedSteps", completedSteps);
    }
    else
    {
        mesh_.metaData().add("completedSteps", completedSteps);
    }
}

bool workflowControls::isStepCompleted() const
{
    const word latestStep = lastCompletedStep();

    if( latestStep.empty() )
        return false;

    const label currVal = workflowSteps_.find(currentStep_)->second;
    const label latestVal = workflowSteps_.find(latestStep)->second;

    if( latestVal == currVal )
        return true;

    return false;
}

bool workflowControls::exitAfterCurrentStep() const
{
    const dictionary& meshDict =
        mesh_.returnTime().lookupObject<dictionary>("meshDict");

    if
    (
        meshDict.found("workflowControls") &&
        meshDict.isDict("workflowControls")
    )
    {
        const dictionary& workflowControls =
            meshDict.subDict("workflowControls");

        if( workflowControls.found("stopAfter") )
        {
            const word exitStep(workflowControls.lookup("stopAfter"));

            if( exitStep == currentStep_ )
                return true;
        }
    }

    return false;
}

word workflowControls::lastCompletedStep() const
{
    if( mesh_.metaData().found("lastStep") )
    {
        const word latestStep(mesh_.metaData().lookup("lastStep"));

        return latestStep;
    }

    return word();
}

DynList<word> workflowControls::completedSteps() const
{
    DynList<word> completedSteps;

    if( mesh_.metaData().found("completedSteps") )
        completedSteps = wordList(mesh_.metaData().lookup("completedSteps"));

    return completedSteps;
}

void workflowControls::clearCompletedSteps()
{
    mesh_.metaData().remove("completedSteps");
    mesh_.metaData().remove("lastStep");
}

bool workflowControls::stopAfterCurrentStep() const
{
    setStepCompleted();

    if( exitAfterCurrentStep() )
    {
        bool writeSuccess(true);

        try
        {
            Info << "Saving mesh generated after step " << currentStep_ << endl;
            mesh_.write();
        }
        catch(...)
        {
            writeSuccess = false;
        }

        returnReduce(writeSuccess, minOp<bool>());

        if( !writeSuccess )
            FatalErrorIn
            (
                "bool workflowControls::stopAfterCurrentStep() const"
            ) << "Mesh was not written on disk" << exit(FatalError);


        std::string message("Stopping after step ");
        message += currentStep_;

        throw message;

        return true;
    }

    return false;
}

bool workflowControls::runAfterCurrentStep() const
{
    if( currentStep_ == restartAfterStep_ )
    {
        try
        {
            Info << "Reading mesh generated after step "
                 << currentStep_ << endl;

            mesh_.read();

            isRestarted_ = true;

            return true;
        }
        catch(...)
        {
            FatalErrorIn
            (
                "bool workflowControls::restartAfterCurrentStep() const"
            ) << "Mesh cannot be loaded. Exitting..." << exit(FatalError);
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

std::map<word, label> workflowControls::populateWorkflowSteps()
{
    std::map<word, label> workflowSteps;
    workflowSteps.insert(std::make_pair(word("start"), 0));
    workflowSteps.insert(std::make_pair(word("templateGeneration"), 1));
    workflowSteps.insert(std::make_pair(word("surfaceTopology"), 2));
    workflowSteps.insert(std::make_pair(word("surfaceProjection"), 4));
    workflowSteps.insert(std::make_pair(word("patchAssignment"), 8));
    workflowSteps.insert(std::make_pair(word("edgeExtraction"), 16));
    workflowSteps.insert(std::make_pair(word("meshOptimisation"), 32));
    workflowSteps.insert(std::make_pair(word("boundaryLayerGeneration"), 64));
    workflowSteps.insert(std::make_pair(word("boundaryLayerRefinement"), 128));

    return workflowSteps;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

workflowControls::workflowControls(polyMeshGen& mesh)
:
    mesh_(mesh),
    currentStep_("start"),
    restartAfterStep_(),
    completedStepsBeforeRestart_(),
    isRestarted_(false)
{
    if( restartRequested() )
    {
        restartAfterStep_ = lastCompletedStep();
        completedStepsBeforeRestart_ = completedSteps();
    }
    else
    {
        clearCompletedSteps();

    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

workflowControls::~workflowControls()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool workflowControls::runCurrentStep(const word& stepName)
{
    if
    (
        completedStepsBeforeRestart_.size() &&
        completedStepsBeforeRestart_.contains(currentStep_) &&
        restartRequested() &&
        !isRestarted_
    )
    {
        Info << "Step " << currentStep_ << " has already been executed" << endl;

        const bool retVal = runAfterCurrentStep();

        //- this step has already been executed
        setStepCompleted();
        currentStep_ = stepName;

        return retVal;
    }
    else if( stopAfterCurrentStep() )
    {
        //- the process shall exit within the stopAfterCurrentStep function
        return false;
    }

    //- check if the requested step exists in the database of steps
    std::map<word, label>::const_iterator it = workflowSteps_.find(stepName);
    if( it == workflowSteps_.end() )
    {
        DynList<word> toc;
        for(it=workflowSteps_.begin();it!=workflowSteps_.end();++it)
            toc.append(it->first);

        FatalErrorIn
        (
            "void workflowControls::setCurrentStep(const word&)"
        ) << "Step " << stepName << " is not a valid name."
          << " Valid step names are " << toc << exit(FatalError);
    }

    setStepCompleted();
    currentStep_ = stepName;

    return true;
}

void workflowControls::workflowCompleted()
{
    if( mesh_.metaData().found("lastStep") )
        mesh_.metaData().remove("lastStep");

    if( mesh_.metaData().found("completedSteps") )
        mesh_.metaData().remove("completedSteps");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
