/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "leastSquaresVectors.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::leastSquaresVectors(const fvMesh& mesh)
:
    MeshObject<fvMesh, leastSquaresVectors>(mesh),
    pVectorsPtr_(NULL),
    nVectorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::~leastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "leastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    pVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();
    const volVectorField& C = mesh().C();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh().nCells(), symmTensor::zero);

    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector d = C[nei] - C[own];
        symmTensor wdd = (1.0/magSqr(d))*sqr(d);

        dd[own] += wdd;
        dd[nei] += wdd;
    }

    forAll(lsP.boundaryField(), patchI)
    {
        const fvPatch& p = mesh().boundary()[patchI];
        const unallocLabelList& faceCells = p.patch().faceCells();

        // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
        const vectorField pd = p.delta();

        forAll(pd, patchFaceI)
        {
            const vector& d = pd[patchFaceI];
            dd[faceCells[patchFaceI]] += (1.0/magSqr(d))*sqr(d);
        }
    }

    // For easy access of neighbour coupled patch field needed for
    // lsN vectors on implicitly coupled boundaries. VV, 18/June/2014
    volSymmTensorField volInvDd
    (
        IOobject
        (
            "volInvDd",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero),
        "zeroGradient"
    );
    symmTensorField& invDd = volInvDd.internalField();
//    invDd = inv(dd);
    invDd = hinv(dd);

    // Evaluate coupled to exchange coupled neighbour field data
    // across coupled boundaries.  HJ, 18/Mar/2015
    volInvDd.boundaryField().evaluateCoupled();

    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector d = C[nei] - C[own];
        scalar magSfByMagSqrd = 1.0/magSqr(d);

        lsP[faceI] = magSfByMagSqrd*(invDd[own] & d);
        lsN[faceI] = -magSfByMagSqrd*(invDd[nei] & d);
    }

    forAll(lsP.boundaryField(), patchI)
    {
        fvsPatchVectorField& patchLsP = lsP.boundaryField()[patchI];
        fvsPatchVectorField& patchLsN = lsN.boundaryField()[patchI];
        const fvPatch& p = mesh().boundary()[patchI];
        const unallocLabelList& faceCells = p.faceCells();

        // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
        const vectorField pd = p.delta();

        if (p.coupled())
        {
            const symmTensorField invDdNei =
                volInvDd.boundaryField()[patchI].patchNeighbourField();

            forAll(pd, patchFaceI)
            {
                const vector& d = pd[patchFaceI];

                patchLsP[patchFaceI] = (1.0/magSqr(d))
                   *(invDd[faceCells[patchFaceI]] & d);

                patchLsN[patchFaceI] = - (1.0/magSqr(d))
                   *(invDdNei[patchFaceI] & d);
            }
        }
        else
        {
            forAll(pd, patchFaceI)
            {
                const vector& d = pd[patchFaceI];

                patchLsP[patchFaceI] = (1.0/magSqr(d))
                   *(invDd[faceCells[patchFaceI]] & d);
            }
        }
    }

    if (debug)
    {
        Info<< "leastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


bool Foam::leastSquaresVectors::movePoints() const
{
    if (debug)
    {
        InfoIn("bool leastSquaresVectors::movePoints() const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}


bool Foam::leastSquaresVectors::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn
        (
            "bool leastSquaresVectors::updateMesh(const mapPolyMesh&) const"
        )   << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}

// ************************************************************************* //
