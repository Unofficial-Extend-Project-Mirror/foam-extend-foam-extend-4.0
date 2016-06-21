/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Class
    numericFlux

Description
    Generic Godunov flux class for coarse level multigrid

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#include "FASNumericFlux.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux>
Foam::FASNumericFlux<Flux>::FASNumericFlux
(
    mgMeshLevel const& meshLevel,
    mgFieldLevel const& fieldLevel
)
:
    meshLevel_(meshLevel),
    fieldLevel_(fieldLevel),
    p_(fieldLevel.pVar()),
    U_(fieldLevel.UVar()),
    T_(fieldLevel.TVar()),
    fineRhoFlux_(fieldLevel.rhoFlux()),
    fineRhoUFlux_(fieldLevel.rhoUFlux()),
    fineRhoEFlux_(fieldLevel.rhoEFlux()),
    rhoFlux_(scalarField(meshLevel.nInternalFaces())),
    rhoUFlux_(vectorField(meshLevel.nInternalFaces())),
    rhoEFlux_(scalarField(meshLevel.nInternalFaces())),
    rhoResidual_(scalarField(meshLevel.nCells())),
    rhoUResidual_(vectorField(meshLevel.nCells())),
    rhoEResidual_(scalarField(meshLevel.nCells()))
{
    rhoFlux_ = 0;
    rhoUFlux_ = vector::zero;
    rhoEFlux_ = 0;
    rhoResidual_ = 0;
    rhoUResidual_ = vector::zero;
    rhoEResidual_ = 0;
    Info<< "mgLevel = " << fieldLevel.level() << endl;
    Info<< "p_ size = " << p_.size() << endl;
    Info<< "U_ size = " << U_.size() << endl;
    Info<< "T_ size = " << T_.size() << endl;
    Info<< "fineRhoFlux_ size = " << fineRhoFlux_.size() << endl;
    Info<< "fineRhoUFlux_ size = " << fineRhoUFlux_.size() << endl;
    Info<< "fineRhoEFlux_ size = " << fineRhoEFlux_.size() << endl;
    Info<< "rhoFlux_ size = " << rhoFlux_.size() << endl;
    Info<< "rhoUFlux_ size = " << rhoUFlux_.size() << endl;
    Info<< "rhoEFlux_ size = " << rhoEFlux_.size() << endl;
    Info<< "rhoResidual_ size = " << rhoResidual_.size() << endl;
    Info<< "rhoUResidual_ size = " << rhoUResidual_.size() << endl;
    Info<< "rhoEResidual_ size = " << rhoEResidual_.size() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux>
inline
void Foam::FASNumericFlux<Flux>::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner     = meshLevel_.owner();
    const unallocLabelList& neighbour = meshLevel_.neighbour();

    // Get the face area vector
    const vectorField& Sf    = meshLevel_.faceAreas();
    const scalarField& magSf = meshLevel_.magFaceAreas();

    const scalarField& Cv = fieldLevel_.CvVar();
    const scalarField& R  = fieldLevel_.R();

    const vectorField& cellCenter = meshLevel_.cellCentres();
    const vectorField& faceCenter = meshLevel_.faceCentres();

    // Calculate fluxes at internal faces
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector deltaRLeft  = faceCenter[faceI] - cellCenter[own];
        vector deltaRRight = faceCenter[faceI] - cellCenter[nei];

        // calculate fluxes with reconstructed primitive variables at faces
        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            p_[own],
            p_[nei],
            U_[own],
            U_[nei],
            T_[own],
            T_[nei],
            R [own],
            R [nei],
            Cv[own],
            Cv[nei],
            Sf[faceI],
            magSf[faceI]
        );
        rhoResidual_[own] += rhoFlux_[faceI];
        rhoResidual_[nei] -= rhoFlux_[faceI];
        rhoUResidual_[own] += rhoUFlux_[faceI];
        rhoUResidual_[nei] -= rhoUFlux_[faceI];
        rhoEResidual_[own] += rhoEFlux_[faceI];
        rhoEResidual_[nei] -= rhoEFlux_[faceI];
    }
    Info<< "Done with interior on level " << fieldLevel_.level() << endl;

    // Update boundary field and values
    forAll(fineRhoFlux_.boundaryField(), patchi)
    {
        unallocLabelList const& owner = meshLevel_.faceCells(patchi);

        const scalarField& pp = fieldLevel_.patchP(patchi);
        const vectorField& pU = fieldLevel_.patchU(patchi);
        const scalarField& pT = fieldLevel_.patchT(patchi);

        const scalarField& pCv = fieldLevel_.patchCv(patchi);
        const scalarField& pR  = fieldLevel_.patchR(patchi);

        const vectorField& pSf    = meshLevel_.patchFaceAreas(patchi);
        const scalarField& pMagSf = meshLevel_.magPatchFaceAreas(patchi);

        if (fieldLevel_.p().boundaryField()[patchi].coupled())
        {
            const scalarField ppLeft  =
                fieldLevel_.p().boundaryField()[patchi].patchInternalField();
            const scalarField ppRight =
                fieldLevel_.p().boundaryField()[patchi].patchNeighbourField();
            const vectorField pULeft  =
                fieldLevel_.U().boundaryField()[patchi].patchInternalField();
            const vectorField pURight =
                fieldLevel_.U().boundaryField()[patchi].patchNeighbourField();
            const scalarField pTLeft  =
                fieldLevel_.T().boundaryField()[patchi].patchInternalField();
            const scalarField pTRight =
                fieldLevel_.T().boundaryField()[patchi].patchNeighbourField();

            forAll(owner, facei)
            {
                label own = owner[facei];
                Flux::evaluateFlux
                (
                    rhoFlux_[facei],
                    rhoUFlux_[facei],
                    rhoEFlux_[facei],
                    ppLeft[facei],
                    ppRight[facei],
                    pULeft[facei],
                    pURight[facei],
                    pTLeft[facei],
                    pTRight[facei],
                    pR[facei],
                    pR[facei],
                    pCv[facei],
                    pCv[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
                rhoResidual_[own]  += rhoFlux_[facei];
                rhoUResidual_[own] += rhoUFlux_[facei];
                rhoEResidual_[own] += rhoEFlux_[facei];
            }
        }
        else
        {
            forAll(owner, facei)
            {
                label own = owner[facei];
                Flux::evaluateFlux
                (
                    rhoFlux_[facei],
                    rhoUFlux_[facei],
                    rhoEFlux_[facei],
                    pp[facei],
                    pp[facei],
                    pU[facei],
                    pU[facei],
                    pT[facei],
                    pT[facei],
                    pR[facei],
                    pR[facei],
                    pCv[facei],
                    pCv[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
                rhoResidual_[own]  += rhoFlux_[facei];
                rhoUResidual_[own] += rhoUFlux_[facei];
                rhoEResidual_[own] += rhoEFlux_[facei];
            }
        }
    }
    Info<< "Done with the boundary on level " << fieldLevel_.level() << endl;
}


// ************************************************************************* //

