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

#include "magLongDelta.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(magLongDelta, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::magLongDelta::magLongDelta(const fvMesh& mesh)
:
    MeshObject<fvMesh, magLongDelta>(mesh),
    magLongDeltaPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::magLongDelta::~magLongDelta()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::magLongDelta::clearOut() const
{
    deleteDemandDrivenData(magLongDeltaPtr_);
}


void Foam::magLongDelta::makeMagLongDistance() const
{
    if (magLongDeltaPtr_)
    {
        FatalErrorIn("void magLongDelta::makeMagLongDistance() const")
            << "Long cell distances already calculated"
            << abort(FatalError);
    }

//     if (debug)
    {
        InfoIn("magLongDelta::makeMagLongDistance()")
            << "Constructing magnitude of long cell distance"
            << endl;
    }

    magLongDeltaPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "magLongDelta",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    );
    surfaceScalarField& mldp = *magLongDeltaPtr_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& Cf = mesh().faceCentres();
    const vectorField& C = mesh().cellCentres();
    const vectorField& Sf = mesh().faceAreas();
    const scalarField& magSf = mesh().magSf();

    forAll (owner, facei)
    {
        // This must be the same as in surfaceInterpolation.C
        scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[owner[facei]]));
        scalar SfdNei = mag(Sf[facei] & (C[neighbour[facei]] - Cf[facei]));
        mldp[facei] = (SfdOwn + SfdNei)/magSf[facei];
    }

    forAll (mldp.boundaryField(), patchi)
    {
        mldp.boundaryField()[patchi] = calcMagLongDistance(patchi);
    }

    if (debug)
    {
        InfoIn("magLongDelta::makeMagLongDistance()")
            << "Finished magnitude of long cell distance"
            << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::magLongDelta::calcMagLongDistance
(
    const label patchi
) const
{
    const fvPatch& p = mesh().boundary()[patchi];

    vectorField d = p.fvPatch::delta();

    if (p.coupled())
    {
        return (mag(p.Sf() & d) + mag(p.Sf() & (p.delta() - d)))/p.magSf();
    }
    else
    {
        return mag(p.Sf() & d)/p.magSf();
    }
}


const Foam::surfaceScalarField& Foam::magLongDelta::magDelta() const
{
    if (!magLongDeltaPtr_)
    {
        makeMagLongDistance();
    }

    return *magLongDeltaPtr_;
}


const Foam::scalarField& Foam::magLongDelta::magDelta
(
    const label patchi
) const
{
    return magDelta().boundaryField()[patchi];
}


bool Foam::magLongDelta::movePoints() const
{
    if (debug)
    {
        InfoIn("bool magLongDelta::movePoints() const")
            << "Clearing long cell distance data" << endl;
    }

    clearOut();

    return true;
}


bool Foam::magLongDelta::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn("bool magLongDelta::updateMesh(const mapPolyMesh&) const")
            << "Clearing long cell distance data" << endl;
    }

    clearOut();

    return true;
}


// ************************************************************************* //
