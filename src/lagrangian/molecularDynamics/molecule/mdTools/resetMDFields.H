if (runTime.outputTime())
{
    allSpeciesN_RU = List< scalarField >
    (
        molecules.pairPotentials().nIds(),
        scalarField (mesh.nCells(), 0.0)
    );

    allSpeciesM_RU = List< scalarField >
    (
        molecules.pairPotentials().nIds(),
        scalarField (mesh.nCells(), 0.0)
    );

    allSpeciesVelocitySum_RU = List< vectorField >
    (
        molecules.pairPotentials().nIds(),
        vectorField (mesh.nCells(), vector::zero)
    );

    allSpeciesVelocityMagSquaredSum_RU = List< scalarField >
    (
        molecules.pairPotentials().nIds(),
        scalarField (mesh.nCells(), 0.0)
    );
}
