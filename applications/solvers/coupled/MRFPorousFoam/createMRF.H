    MRFZones mrfZones(mesh);
    mrfZones.correctBoundaryVelocity(U);

    bool addMRF = false;

    if (!mrfZones.empty())
    {
        addMRF = true;
    }
