    if (runTime.outputTime())
    {
        OFstream nAveragingStepsFile
        (
            runTime.path()/runTime.timeName()/"nAveragingSteps.raw"
        );

        if (nAveragingStepsFile.good())
        {
            nAveragingStepsFile << nAveragingSteps << endl;
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Cannot open file "
                << runTime.path()/runTime.timeName()/"nAveragingSteps.raw"
                << abort(FatalError);
        }
    }
