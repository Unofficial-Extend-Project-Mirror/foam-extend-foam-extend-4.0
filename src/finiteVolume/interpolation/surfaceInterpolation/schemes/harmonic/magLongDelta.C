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
    magLongDeltaPtr_(NULL),
    magLongDeltaBnd_(mesh.boundary().size(), NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::magLongDelta::~magLongDelta()
{
    clearData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::magLongDelta::clearData() const
{
    if (magLongDeltaPtr_)
    {
        delete magLongDeltaPtr_;
    }
    else
    {
        forAll (magLongDeltaBnd_, i)
        {
            if (magLongDeltaBnd_[i])
            {
                delete magLongDeltaBnd_[i];
            }
        }
    }
}


void Foam::magLongDelta::makeMagLongDistance() const
{
    if (debug)
    {
        Info<< "magLongDelta::makeMagLongDistance() :"
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
        const fvPatch& p = mesh().boundary()[patchi];

        if (p.coupled())
        {
            if (!magLongDeltaBnd_[patchi])
            {
                makeMagLongDistance(patchi);
            }

            mldp.boundaryField()[patchi] = *magLongDeltaBnd_[patchi];
            delete magLongDeltaBnd_[patchi];
            magLongDeltaBnd_[patchi] = &mldp.boundaryField()[patchi];
        }
        else
        {
            delete magLongDeltaBnd_[patchi];
        }
    }

    if (debug)
    {
        Info<< "magLongDelta::makeMagLongDistance() :"
            << "Finished magnitude of long cell distance"
            << endl;
    }
}


void Foam::magLongDelta::makeMagLongDistance(label patchi) const
{
    if (debug)
    {
        Info<< "magLongDelta::makeMagLongDistanceBnd(label patchi) :"
            << "Constructing magnitude of long cell distance"
            << endl;
    }

    const fvPatch& p = mesh().boundary()[patchi];

    vectorField d = p.fvPatch::delta();

    magLongDeltaBnd_[patchi] =
        new scalarField
        (
            (mag(p.Sf() & d) + mag(p.Sf() & (p.delta() - d)))/p.magSf()
        );

    if (debug)
    {
        Info<< "magLongDelta::makeMagLongDistanceBnd(label patchi) :"
            << "Finished magnitude of long cell distance"
            << endl;
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
    if (!mesh().boundary()[patchi].coupled())
    {
        FatalErrorIn
        (
            "const Foam::scalarField& Foam::magLongDelta::magDelta("
            "const label patchi) const"
        )   << "patch is not a coupled. Cannot calculate long distance"
            << abort(FatalError);
    }

    if (!magLongDeltaBnd_[patchi])
    {
        makeMagLongDistance(patchi);
    }

    return *magLongDeltaBnd_[patchi];
}


bool Foam::magLongDelta::movePoints() const
{
    if (debug)
    {
        InfoIn("bool magLongDelta::movePoints() const")
            << "Clearing long cell distance data" << endl;
    }

    clearData();

    magLongDeltaBnd_.setSize(mesh().boundary().size(), NULL);

    return true;
}


bool Foam::magLongDelta::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn("bool magLongDelta::updateMesh(const mapPolyMesh&) const")
            << "Clearing long cell distance data" << endl;
    }

    clearData();

    magLongDeltaBnd_.setSize(mesh().boundary().size(), NULL);

    return true;
}


// ************************************************************************* //
