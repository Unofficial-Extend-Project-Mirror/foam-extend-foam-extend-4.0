{
    label sizeByTwo =  mesh.boundary()[cohesivePatchID].size()/2;
    vectorField Cczm = mesh.boundary()[cohesivePatchID].Cf();

    for(label i = 0; i < sizeByTwo; i++)
    {
        Cczm[i] -= Cczm[sizeByTwo + i];
        Cczm[sizeByTwo + i] = Cczm[i];
    }

    scalar error = sum(mag(Cczm));

    if(error < SMALL)
    {
        Info << "Cohesive patch face ordering: OK" << endl;
    }
    else
    {
        Info << "Cohesive patch face ordering: error (" << error << ")"
            << endl;
    }
}
