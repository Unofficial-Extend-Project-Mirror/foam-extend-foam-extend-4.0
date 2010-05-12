    mkDir(runTime.path()/"positions"/runTime.timeName());

    OFstream positionsFile
    (
        runTime.path()
        /"positions"
        /runTime.timeName()
        /"leftAndRightPosition"
    );

    positionsFile << "Time" << tab;

    positionsFile << "Left" << tab << "Right" << endl;


    // Find left and right point label in the patch localPoints field

    const pointField& locPoints =
        mesh.boundaryMesh()[interface.aPatchID()].localPoints();

    label leftPointLabel = -1;
    label rightPointLabel = -1;

    forAll(locPoints, pointI)
    {
        if(leftPointLabel == -1 || rightPointLabel == -1)
        {
            if(mag(locPoints[pointI].x()) < SMALL)
            {
                leftPointLabel = pointI;
            }
            else if(mag(locPoints[pointI].x() - 1.0) < SMALL)
            {
                rightPointLabel = pointI;
            }
        }
    }

    if(leftPointLabel == -1 || rightPointLabel == -1)
    {
        Info << "Can't find left and right points labels!";
        return 0;
    }

// scalar dissipation = 0.0;
