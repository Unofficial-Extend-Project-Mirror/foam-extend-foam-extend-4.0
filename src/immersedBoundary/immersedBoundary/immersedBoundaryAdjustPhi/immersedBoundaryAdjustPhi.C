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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryAdjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "immersedBoundaryFvPatchFields.H"
#include "immersedBoundaryFvsPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::immersedBoundaryAdjustPhi
(
    surfaceScalarField& phi,
    volVectorField& U
)
{
    const fvMesh& mesh = phi.mesh();

    // If the mesh is moving, adjustment needs to be calculated on
    // relative fluxes.  HJ, 13/Feb/2009
    if (mesh.moving())
    {
        fvc::makeRelative(phi, U);
    }

    forAll (phi.boundaryField(), patchI)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchI];

        if (isA<immersedBoundaryFvPatchVectorField>(Up))
        {
            if (Up.fixesValue())
            {
                // Found immersed boundary path which fixes value.
                // Correction is necessary

                // Cast into immersedBoundary types
                const immersedBoundaryFvPatchVectorField& ibU =
                    refCast<const immersedBoundaryFvPatchVectorField>(Up);

                const immersedBoundaryFvPatch& ibPatch = ibU.ibPatch();

                // Complete immersed boundary triangular surface is present
                // on all processors: no need for parallel reduction

                // Sum the flux through the immersed boundary
                const scalar triFlux = sum(ibU.refValue() & ibPatch.triSf());

                scalarField& phiInternal = phi.internalField();

                scalar fluxIn = 0;
                scalar fluxOut = 0;
                scalar fixedFlux = 0;

                // Get all IB faces
                const labelList& ibFaces = ibPatch.ibFaces();
                const boolList& ibFaceFlips = ibPatch.ibFaceFlips();

                forAll (ibFaces, faceI)
                {
                    const label curFace = ibFaces[faceI];
                    const bool curFlip = ibFaceFlips[faceI];

                    if (mesh.isInternalFace(curFace))
                    {
                        const scalar curFlux = phiInternal[curFace];

                        if (!curFlip)
                        {
                            // Face points out of the live cell
                            if (curFlux >= 0)
                            {
                                // Flux out of the live cell
                                fluxOut += curFlux;
                            }
                            else
                            {
                                // Flux into the live cell
                                fluxIn -= curFlux;
                            }
                        }
                        else
                        {
                            // Face points into the live cell: flip it
                            if (curFlux >= 0)
                            {
                                // Flux into of the live cell
                                fluxIn += curFlux;
                            }
                            else
                            {
                                // Flux out the live cell
                                fluxOut -= curFlux;
                            }
                        }
                    }
                    else
                    {
                        const label patchID =
                            mesh.boundaryMesh().whichPatch(curFace);
                        const label faceID =
                            mesh.boundaryMesh()[patchID].whichFace(curFace);

                        const scalar curFlux =
                            phi.boundaryField()[patchID][faceID];

                        // Note: only coupled patches may carry flux
                        // In order to avoid double summation and
                        // inconsistencies in the correction,
                        // coupled face fluxes will NOT be corrected,
                        // but only accounted for in the summation.
                        if (mesh.boundaryMesh()[patchID].coupled())
                        {
                            // Only do the master side; slave will
                            // be handled on the other side of the couple
                            if (!curFlip)
                            {
                                fixedFlux += curFlux;
                            }
                        }
                    }
                }

                reduce(fluxIn, sumOp<scalar>());
                reduce(fluxOut, sumOp<scalar>());
                reduce(fixedFlux, sumOp<scalar>());

                scalar imbalance = (fluxIn - fluxOut + fixedFlux) - triFlux;

                if (fvMesh::debug)
                {
                    Info<< "triFlux = " << triFlux
                        << " fluxIn = " << fluxIn << " fluxOut = " << fluxOut
                        << " fixedFlux = " << fixedFlux
                        << " imbalance = " << imbalance
                        << endl;
                }

                scalar massCorr = 1.0;

                if (mag(imbalance) > SMALL)
                {
                    // Scaling required: scale to match the smaller of two
                    // fluxes
                    if (fluxIn > fluxOut)
                    {
                        // Scale down incoming flux
                        // Note change of sign: imbalance is negative
                        massCorr = 1 - imbalance/(fluxIn + SMALL);

                        if (fvMesh::debug)
                        {
                            Info<< "Scaling down incoming flux with factor = "
                                << massCorr << endl;
                        }

                        scalar newFluxIn = 0;

                        // Visit all incoming flux faces and re-scale the flux
                        forAll (ibFaces, faceI)
                        {
                            const label curFace = ibFaces[faceI];
                            const bool curFlip = ibFaceFlips[faceI];

                            if (mesh.isInternalFace(curFace))
                            {
                                // Take reference to current flux
                                scalar& curFlux = phiInternal[curFace];

                                if (!curFlip)
                                {
                                    // Face points out of the live cell
                                    if (curFlux < 0)
                                    {
                                        // Flux out of the live cell
                                        curFlux *= massCorr;
                                        newFluxIn -= curFlux;
                                    }
                                }
                                else
                                {
                                    // Face points into the live cell: flip it
                                    if (curFlux >= 0)
                                    {
                                        // Flux out the live cell
                                        curFlux *= massCorr;
                                        newFluxIn += curFlux;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        // Scale down outgoing flux
                        massCorr = 1 + imbalance/(fluxOut + SMALL);

                        if (fvMesh::debug)
                        {
                            Info<< "Scaling down outgoing flux with factor = "
                                << massCorr << endl;
                        }

                        scalar newFluxOut = 0;

                        // Visit all outgoing flux faces and re-scale the flux
                        forAll (ibFaces, faceI)
                        {
                            const label curFace = ibFaces[faceI];
                            const bool curFlip = ibFaceFlips[faceI];

                            if (mesh.isInternalFace(curFace))
                            {
                                // Take reference to current flux
                                scalar& curFlux = phiInternal[curFace];

                                if (!curFlip)
                                {
                                    // Face points out of the live cell
                                    if (curFlux >= 0)
                                    {
                                        // Flux out of the live cell
                                        curFlux *= massCorr;
                                        newFluxOut += curFlux;
                                    }
                                }
                                else
                                {
                                    // Face points into the live cell: flip it
                                    if (curFlux < 0)
                                    {
                                        // Flux out the live cell
                                        curFlux *= massCorr;
                                        newFluxOut -= curFlux;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    // If the mesh is moving, adjustment needs to be calculated on
    // relative fluxes.  Now reverting to absolute fluxes.  HJ, 13/Feb/2009
    if (mesh.moving())
    {
        fvc::makeAbsolute(phi, U);
    }
}


// ************************************************************************* //
