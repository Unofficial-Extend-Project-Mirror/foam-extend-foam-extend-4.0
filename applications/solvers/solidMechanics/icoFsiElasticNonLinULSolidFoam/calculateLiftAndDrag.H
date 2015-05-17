{
    //- calculate lift and drag for plate and cylinder due to pressure (ignore shear stress)
    label plateID = mesh.boundaryMesh().findPatchID("plate");
    label cylinderID = mesh.boundaryMesh().findPatchID("cylinder");

    if(plateID == -1 || cylinderID == -1)
    {
        FatalError << "\n Cannot find the plate patch or the cylinder"
            << " patch to calculate lift and drag!"
            << exit(FatalError);
    }

    scalar lift = 0;
    scalar drag = 0;

    const vectorField& Sfp = mesh.boundary()[plateID].Sf();
    forAll(p.boundaryField()[plateID], facei)
    {
        vector faceForce = p.boundaryField()[plateID][facei] * Sfp[facei];
        lift += vector(0,1,0) & faceForce;
        drag += vector(1,0,0) & faceForce;
    }

    const vectorField& Sfc = mesh.boundary()[cylinderID].Sf();
    forAll(p.boundaryField()[cylinderID], facei)
    {
        vector faceForce = p.boundaryField()[cylinderID][facei] * Sfc[facei];
        lift += vector(0,1,0) & faceForce;
        drag += vector(1,0,0) & faceForce;
    }

    scalar width = 0.050668;

    Info<< "Total lift on the cylinder and plate boundaries is " << lift << " N, per unit width is " << (lift/width) << " N\n"
        << "Total drag on the cylinder and plate boundaries is " << drag << " N, per unit width is " << (drag/width) << " N\n"
        << endl;
}
