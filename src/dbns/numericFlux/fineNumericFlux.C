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

\*---------------------------------------------------------------------------*/

#include "fineNumericFlux.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux, class Limiter>
Foam::fineNumericFlux<Flux, Limiter>::fineNumericFlux
(
    const mgMeshLevel& meshLevel,
    const mgFieldLevel& fieldLevel
    basicThermo& thermo,
)
:
    meshLevel_(meshLevel),
    fieldLevel_(fieldLevel),
    p_(fieldLevel.p()),
    U_(fieldLevel.U()),
    T_(fieldLevel.T()),
    thermo_(thermo),
    rhoFlux_(fieldLevel.rhoFlux()),
    rhoUFlux_(fieldLevel.rhoUFlux()),
    rhoEFlux_(fieldLevel.rhoEFlux())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux, class Limiter>
void Foam::fineNumericFlux<Flux, Limiter>::computeInteriorFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = meshLevel_.owner();
    const unallocLabelList& neighbour = meshLevel_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh().Sf();
    const surfaceScalarField& magSf = mesh().magSf();

    const scalarField& Cv = fieldLevel_.CvVar();
    const scalarField& R  = fieldLevel_.R();

    // Get gradients
    // Coupled patch update on gradients moved into gradScheme.C
    // HJ, 22/Apr/2016;

    // Changed return type for gradient cacheing.  HJ, 22/Apr/2016
    const tmp<volVectorField> tgradP = fvc::grad(p_);
    const volVectorField& gradP = tgradP();

    const tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    const tmp<volVectorField> tgradT = fvc::grad(T_);
    const volVectorField& gradT = tgradT();

    MDLimiter<scalar, Limiter> scalarPLimiter
    (
        this->p_,
        gradP
    );

    MDLimiter<vector, Limiter> vectorULimiter
    (
        this->U_,
        gradU
    );

    MDLimiter<scalar, Limiter> scalarTLimiter
    (
        this->T_,
        gradT
    );

    // Get limiters
    const volScalarField& pLimiter = scalarPLimiter.phiLimiter();
    const volVectorField& ULimiter = vectorULimiter.phiLimiter();
    const volScalarField& TLimiter = scalarTLimiter.phiLimiter();

    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
        const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

        // calculate fluxes with reconstructed primitive variables at faces
        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            p_[own] + pLimiter[own]*(deltaRLeft & gradP[own]),
            p_[nei] + pLimiter[nei]*(deltaRRight & gradP[nei]),
            U_[own] + cmptMultiply(ULimiter[own], (deltaRLeft & gradU[own])),
            U_[nei] + cmptMultiply(ULimiter[nei], (deltaRRight & gradU[nei])),
            T_[own] + TLimiter[own]*(deltaRLeft & gradT[own]),
            T_[nei] + TLimiter[nei]*(deltaRRight & gradT[nei]),
            R[own],
            R[nei],
            Cv[own],
            Cv[nei],
            Sf[faceI],
            magSf[faceI]
        );
    }

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryField()[patchi];

        const scalarField& pp = fieldLevel_.patchP(patchi);
        const vectorField& pU = fieldLevel_.patchU(patchi);
        const scalarField& pT = fieldLevel_.patchT(patchi);

        const scalarField& pCv = fieldLevel_.patchCv(patchi);
        const scalarField& pR = fieldLevel_.patchR(patchi);

        // Gradients
        const fvPatchVectorField& pGradP = gradP.boundaryField()[patchi];
        const fvPatchTensorField& pGradU = gradU.boundaryField()[patchi];
        const fvPatchVectorField& pGradT = gradT.boundaryField()[patchi];

        // Limiters
        const fvPatchScalarField& pPatchLim = pLimiter.boundaryField()[patchi];
        const fvPatchVectorField& UPatchLim = ULimiter.boundaryField()[patchi];
        const fvPatchScalarField& TPatchLim = TLimiter.boundaryField()[patchi];

        // Face areas
        const vectorField& pSf  = meshLevel_.patchFaceAreas(patchi);
        const scalarField& pMagSf = meshLevel_.magPatchFaceAreas(patchi);

        if (p_.boundaryField()[patchi].coupled())
        {
            // Coupled patch
            const scalarField ppLeft  =
                p_.boundaryField()[patchi].patchInternalField();

            const scalarField ppRight =
                p_.boundaryField()[patchi].patchNeighbourField();

            const vectorField pULeft  =
                U_.boundaryField()[patchi].patchInternalField();

            const vectorField pURight =
                U_.boundaryField()[patchi].patchNeighbourField();

            const scalarField pTLeft  =
                T_.boundaryField()[patchi].patchInternalField();

            const scalarField pTRight =
                T_.boundaryField()[patchi].patchNeighbourField();

            // Gradients
            const vectorField pgradPLeft = pGradP.patchInternalField();
            const vectorField pgradPRight = pGradP.patchNeighbourField();

            const tensorField pgradULeft = pGradU.patchInternalField();
            const tensorField pgradURight = pGradU.patchNeighbourField();

            const vectorField pgradTLeft = pGradT.patchInternalField();
            const vectorField pgradTRight = pGradT.patchNeighbourField();

            // Geometry
            vectorField pDeltaRLeft = meshLevel_.patchDeltaR(patchi);
            vectorField pDdeltaRRight =
                meshLevel_.patchDeltaRNeighbour(patchi);

            // Limiters

            const scalarField ppLimiterLeft = pPatchLim.patchInternalField();
            const scalarField ppLimiterRight = pPatchLim.patchNeighbourField();

            const vectorField pULimiterLeft = UPatchLim.patchInternalField();
            const vectorField pULimiterRight = UPatchLim.patchNeighbourField();

            const scalarField pTLimiterLeft = TPatchLim.patchInternalField();
            const scalarField pTLimiterRight = TPatchLim.patchNeighbourField();

            forAll (pp, facei)
            {
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],

                    ppLeft[facei]
                  + ppLimiterLeft[facei]*
                    (pDeltaRLeft[facei] & pgradPLeft[facei]),

                    ppRight[facei]
                  + ppLimiterRight[facei]*
                    (pDdeltaRRight[facei] & pgradPRight[facei]),

                    pULeft[facei]
                  + cmptMultiply
                    (
                        pULimiterLeft[facei],
                        pDeltaRLeft[facei] & pgradULeft[facei]
                    ),

                    pURight[facei]
                  + cmptMultiply
                    (
                        pULimiterRight[facei],
                        pDdeltaRRight[facei] & pgradURight[facei]
                    ),

                    pTLeft[facei]
                  + pTLimiterLeft[facei]*
                    (pDeltaRLeft[facei] & pgradTLeft[facei]),

                    pTRight[facei]
                  + pTLimiterRight[facei]*
                    (pDdeltaRRight[facei] & pgradTRight[facei]),

                    pR[facei],
                    pR[facei],
                    pCv[facei],
                    pCv[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
        else
        {
            forAll (pp, facei)
            {
                // Calculate fluxes
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
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
            }
        }
    }
}


// ************************************************************************* //
