/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Global
    wallViscosity

Description
    Calculate wall viscosity from wall-functions.

\*---------------------------------------------------------------------------*/

{
    scalar Cmu25 = pow(Cmu_.value(), 0.25);

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            const scalarField& rhow = rho_.boundaryField()[patchi];

            const scalarField& muw = mu().boundaryField()[patchi];
            scalarField& mutw = mut_.boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*RASModel::y_[patchi][facei]
                   *sqrt(k_[faceCelli])/(muw[facei]/rhow[facei]);

                if (yPlus > yPlusLam_)
                {
                    mutw[facei] =
                        muw[facei]
                       *(yPlus*kappa_.value()/log(E_.value()*yPlus) - 1);
                }
                else
                {
                    mutw[facei] = 0.0;
                }
            }
        }
    }
}


// ************************************************************************* //
