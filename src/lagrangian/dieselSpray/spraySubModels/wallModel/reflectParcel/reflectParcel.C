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

#include "reflectParcel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reflectParcel, 0);

addToRunTimeSelectionTable
(
    wallModel,
    reflectParcel,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
reflectParcel::reflectParcel
(
    const dictionary& dict,
    const volVectorField& U,
    spray& sm
)
:
    wallModel(dict, U, sm),
    U_(U),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    elasticity_(readScalar(coeffsDict_.lookup("elasticity")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reflectParcel::~reflectParcel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return 'keepParcel'
bool reflectParcel::wallTreatment
(
    parcel& p,
    const label globalFacei
) const
{
    label patchi = p.patch(globalFacei);
    label facei = p.patchFace(patchi, globalFacei);

    const polyMesh& mesh = spray_.mesh();

    if (mesh_.boundaryMesh()[patchi].isWall())
    {
        // wallNormal defined to point outwards of domain
        vector Sf = mesh_.Sf().boundaryField()[patchi][facei];
        Sf /= mag(Sf);

        if (!mesh.moving())
        {
            // static mesh
            scalar Un = p.U() & Sf;

            if (Un > 0)
            {
                p.U() -= (1.0 + elasticity_)*Un*Sf;
            }

        }
        else
        {
            // moving mesh
            vector Ub1 = U_.boundaryField()[patchi][facei];
            vector Ub0 = U_.oldTime().boundaryField()[patchi][facei];

            scalar dt = spray_.runTime().deltaT().value();
            const vectorField& oldPoints = mesh.oldPoints();

            const vector& Cf1 = mesh.faceCentres()[globalFacei];

            vector Cf0 = mesh.faces()[globalFacei].centre(oldPoints);
            vector Cf = Cf0 + p.stepFraction()*(Cf1 - Cf0);
            vector Sf0 = mesh.faces()[globalFacei].normal(oldPoints);

            // for layer addition Sf0 = vector::zero and we use Sf
            if (mag(Sf0) > SMALL)
            {
                Sf0 /= mag(Sf0);
            }
            else
            {
                Sf0 = Sf;
            }

            scalar magSfDiff = mag(Sf - Sf0);

            vector Ub = Ub0 + p.stepFraction()*(Ub1 - Ub0);

            if (magSfDiff > SMALL)
            {
                // rotation + translation
                vector Sfp = Sf0 + p.stepFraction()*(Sf - Sf0);

                vector omega = Sf0 ^ Sf;
                scalar magOmega = mag(omega);
                omega /= magOmega+SMALL;

                scalar phiVel = ::asin(magOmega)/dt;

                scalar dist = (p.position() - Cf) & Sfp;
                vector pos = p.position() - dist*Sfp;
                vector vrot = phiVel*(omega ^ (pos - Cf));

                vector v = Ub + vrot;

                scalar Un = ((p.U() - v) & Sfp);

                if (Un > 0.0)
                {
                    p.U() -= (1.0 + elasticity_)*Un*Sfp;
                }
            }
            else
            {
                // translation
                vector Ur = p.U() - Ub;
                scalar Urn = Ur & Sf;
                /*
                if (mag(Ub-v) > SMALL)
                {
                    Info << "reflectParcel:: v = " << v
                        << ", Ub = " << Ub
                        << ", facei = " << facei
                        << ", patchi = " << patchi
                        << ", globalFacei = " << globalFacei
                        << ", name = " << mesh_.boundaryMesh()[patchi].name()
                        << endl;
                }
                    */
                if (Urn > 0.0)
                {
                    p.U() -= (1.0 + elasticity_)*Urn*Sf;
                }
            }
        }

    }
    else
    {
        FatalError
            << "bool reflectParcel::wallTreatment(parcel& parcel) const "
                << " parcel has hit a boundary "
                << mesh_.boundary()[patchi].type()
                << " which not yet has been implemented."
            << abort(FatalError);
    }
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
