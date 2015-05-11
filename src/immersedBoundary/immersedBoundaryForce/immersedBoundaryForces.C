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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryForces.H"
#include "immersedBoundaryFvPatch.H"
#include "immersedBoundaryFvPatchFields.H"
#include "immersedBoundaryVelocityWallFunctionFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "Time.H"

using namespace Foam::incompressible::RASModels;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryForces, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryForces::immersedBoundaryForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    forces
    (
        name,
        obr,
        dict,
        loadFromFiles
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryForces::~immersedBoundaryForces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::immersedBoundaryForces::forcesMoments
Foam::immersedBoundaryForces::calcForcesMoment() const
{
    forcesMoments fm
    (
        pressureViscous(vector::zero, vector::zero),
        pressureViscous(vector::zero, vector::zero)
    );

    if (directForceDensity_)
    {
        const volVectorField& fD = obr_.lookupObject<volVectorField>(fDName_);

        const fvMesh& mesh = fD.mesh();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            // Check and cast into immersed boundary type
            if
            (
                isA<immersedBoundaryFvPatchVectorField>
                (
                    fD.boundaryField()[patchI]
                )
            )
            {
                // Found immersed boundary patch and field.
                // Cast into immersed boundary type
                const immersedBoundaryFvPatch& ibPatch =
                    refCast<const immersedBoundaryFvPatch>
                    (
                        mesh.boundary()[patchI]
                    );

                const immersedBoundaryFvPatchVectorField& fDpatch =
                    refCast<const immersedBoundaryFvPatchVectorField>
                    (
                        fD.boundaryField()[patchI]
                    );

                // Get face area vectors for triangles
                const vectorField& Sfb = ibPatch.triSf();
                scalarField sA = mag(Sfb);

                // Calculate distance for triangles
                vectorField Md = ibPatch.triCf() - CofR_;

                // Old treatment:
                // Normal force =
                // surfaceNormal*(surfaceUnitNormal & forceDensity)
                // The first operation will be done on ibPoints, the data will
                // then be distributed onto the ib surface
                // for surface integration
//                vectorField fN =
//                    Sfb*
//                    ibPatch.toTriFaces
//                    (
//                        ibPatch.ibNormals() & fDpatch.ibValue()
//                    );

                // New treatment: normal force calculated on triangles
                // Damir Rigler, 30/Apr/2014
                vectorField fN = Sfb/sA*(Sfb & fDpatch.triValue());

                fm.first().first() += sum(fN);
                fm.second().first() += sum(Md ^ fN);

                // Tangential force (total force minus normal fN)
                vectorField fT = sA*fDpatch.triValue() - fN;

                fm.first().second() += sum(fT);
                fm.second().second() += sum(Md ^ fT);
            }
        }
    }
    else
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        const fvMesh& mesh = U.mesh();

        // Stress from turbulence model.  Note: this does not account
        // for wall functions on the Immersed Boundary: use wallShearStress
        // from Immersed Boundary U wall patch instead.
        // HJ, 11/Aug/2014
//         volSymmTensorField stress = devRhoReff();

        // Scale pRef by density for incompressible simulations
        scalar pRef = pRef_/rho(p);

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            // Check and cast into immersed boundary type
            if
            (
                isA<immersedBoundaryFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
            {
                // Found immersed boundary patch and field.
                // Cast into immersed boundary type
                const immersedBoundaryFvPatch& ibPatch =
                    refCast<const immersedBoundaryFvPatch>
                    (
                        mesh.boundary()[patchI]
                    );

                // Get immersed boundary velocity.  Used to access wall
                // shear stress
                const immersedBoundaryVelocityWallFunctionFvPatchVectorField&
                    UPatch =
                    refCast
                    <
                        const
                        immersedBoundaryVelocityWallFunctionFvPatchVectorField
                    >
                    (
                        U.boundaryField()[patchI]
                    );

//                 const immersedBoundaryFvPatchSymmTensorField stressPatch =
//                     refCast<const immersedBoundaryFvPatchSymmTensorField>
//                     (
//                         stress.boundaryField()[patchI]
//                     );

                const immersedBoundaryFvPatchScalarField pPatch =
                    refCast<const immersedBoundaryFvPatchScalarField>
                    (
                        p.boundaryField()[patchI]
                    );

                // Get face area vectors for triangles
                const vectorField& Sfb = ibPatch.triSf();
                scalarField sA = mag(Sfb);

                // Calculate distance for triangles
                vectorField Md = ibPatch.triCf() - CofR_;

                // Pressure force is an integral of interpolated pressure
                // on triangular faces
                vectorField pf = Sfb*(pPatch.triValue() - pRef);

                fm.first().first() += rho(p)*sum(pf);
                fm.second().first() += rho(p)*sum(Md ^ pf);

                // Old treatment:
                // Shear force is obtained from velocity wall functions
                // and integrated on triangular faces
                vectorField vf =
                    sA*ibPatch.toTriFaces(UPatch.wallShearStress());

                // New treatment: normal force calculated on triangles
                // Damir Rigler, 30/Apr/2014
//                 vectorField vfOld = (Sfb & stressPatch.triValue());
//                 Info<< "Force difference " << sumMag(vf - vfOld) << endl;

                fm.first().second() += sum(vf);
                fm.second().second() += sum(Md ^ vf);
            }
        }
    }

    reduce(fm, sumOp());

    return fm;
}


// ************************************************************************* //
