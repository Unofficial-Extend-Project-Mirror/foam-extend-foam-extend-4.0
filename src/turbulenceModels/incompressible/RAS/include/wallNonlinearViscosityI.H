/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Global
    wallNonlinearViscosity

Description
    Calculate wall viscosity for non-linear models

\*---------------------------------------------------------------------------*/

{
    const fvPatchList& patches = mesh_.boundary();

    const scalar yPlusLam = this->yPlusLam(kappa_.value(), E_.value());

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (curPatch.isWall())
        {
            const scalarField& nuw = nu().boundaryField()[patchi];
            scalarField& nutw = nut_.boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                //- Using local Cmu
                scalar Cmu25 = pow(Cmu_[faceCelli], 0.25);

                scalar yPlus =
                    Cmu25*y_[patchi][facei]*sqrt(k_[faceCelli])/nuw[facei];

                if (yPlus > yPlusLam)
                {
                    nutw[facei] =
                        nuw[facei]
                       *(yPlus*kappa_.value()/log(E_.value()*yPlus) - 1);
                }
                else
                {
                    nutw[facei] = 0.0;
                }
            }
        }
    }
}


// ************************************************************************* //
