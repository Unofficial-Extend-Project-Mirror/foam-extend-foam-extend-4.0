{
    scalar Cmu25 = ::pow(Cmu, 0.25);

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isType<wallFvPatch>(currPatch))
        {
            scalarField& nutw = nutb.boundaryField()[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                // calculate yPlus
                scalar yPlus =
                    Cmu25*y[patchi][facei]
                   *::sqrt(k[faceCelli])
                   /nub.value();

                if (yPlus > 11.6)
                {
                    nutw[facei] =
                        yPlus*nub.value()*kappa
                       /::log(E*yPlus)
                      - nub.value();
                }
                else
                {
                    nutw[facei] = 0.0;
                }
            }
        }
    }
}
