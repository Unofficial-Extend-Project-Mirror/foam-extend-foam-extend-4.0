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

#include "adjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorFvsPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "fvc.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Note: all tolerances need to be global statics: closedDomainTol
// taken outside of the function.  HJ, 25/Dec/2015
static const Foam::debug::tolerancesSwitch closedDomainTol
(
    "closedDomainTol",
    1e-10
);

bool Foam::adjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& p
)
{
    if (p.needReference())
    {
        // Removed updateCoeffs.  HJ, 9/May/2016

        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;

        // If the mesh is moving, adjustment needs to be calculated on
        // relative fluxes.  HJ, 13/Feb/2009
        if (phi.mesh().moving())
        {
            fvc::makeRelative(phi, U);
        }

        forAll (phi.boundaryField(), patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = phi.boundaryField()[patchi];

            // Bug fix: All coupled patches should also be unaffected
            // HJ, 12/Feb/2010
            if (!Up.coupled())
            {
                if
                (
                    Up.fixesValue()
                 && !isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    forAll (phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }
                }
                else
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            adjustableMassOut += phip[i];
                        }
                    }
                }
            }
        }


        //HR 16.03.10: Bug fix for moving meshes. Changing domain volume needs
        // to be taken into account.
        if (phi.mesh().moving())
        {
            dimensionedScalar Vdiff =
                sum(phi.mesh().V()) - sum(phi.mesh().V0());

            fixedMassOut += Vdiff.value()/phi.time().deltaT().value();
        }

        reduce(massIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());

        scalar massCorr = 1.0;

        if (mag(adjustableMassOut) > SMALL)
        {
            massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }
        else if
        (
            mag(fixedMassOut - massIn)
          > closedDomainTol()*Foam::max(1.0, mag(massIn))
        )
        {
            phi.write();
            U.write();
            p.write();

            // Cannot adjust
            FatalErrorIn
            (
                "adjustPhi\n"
                "(\n"
                "    surfaceScalarField& phi,\n"
                "    const volVectorField& U,\n"
                "    const volScalarField& p\n"
                ")"
            )   << "Continuity error cannot be removed by adjusting the"
                " outflow.\nPlease check the velocity boundary conditions"
                " and/or run potentialFoam to initialise the outflow." << nl
                << "Specified mass inflow   : " << massIn << nl
                << "Specified mass outflow  : " << fixedMassOut << nl
                << "Difference              : "
                << mag(fixedMassOut - massIn) << nl
                << "Adjustable mass outflow : " << adjustableMassOut << nl
                << exit(FatalError);
        }

        if (fvMesh::debug)
        {
            Info<< "bool Foam::adjustPhi(...) massIn: " << massIn
                << " fixedMassOut: " << fixedMassOut
                << " adjustableMassOut: " << adjustableMassOut
                << " mass corr: " << massCorr
                << endl;
        }

        forAll (phi.boundaryField(), patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = phi.boundaryField()[patchi];

            // Bug fix: All coupled patches should also be unaffected
            // HJ, 12/Feb/2010
            if (!Up.coupled())
            {
                if
                (
                    !Up.fixesValue()
                 || isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    forAll (phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                    }
                }
            }
        }

        // If the mesh is moving, adjustment needs to be calculated on
        // relative fluxes.  Now reverting to absolute fluxes.  HJ, 13/Feb/2009
        if (phi.mesh().moving())
        {
            fvc::makeAbsolute(phi, U);
        }

        return mag(massIn) < SMALL
            && mag(fixedMassOut) < SMALL
            && mag(adjustableMassOut) < SMALL;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
