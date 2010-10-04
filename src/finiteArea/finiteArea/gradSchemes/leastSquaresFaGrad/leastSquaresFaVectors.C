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

#include "leastSquaresFaVectors.H"
#include "edgeFields.H"
#include "areaFields.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresFaVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresFaVectors::leastSquaresFaVectors(const faMesh& mesh)
:
    MeshObject<faMesh, leastSquaresFaVectors>(mesh),
    pVectorsPtr_(NULL),
    nVectorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresFaVectors::~leastSquaresFaVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresFaVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "leastSquaresFaVectors::makeLeastSquaresVectors() :"
            << "Constructing finite area least square gradient vectors"
            << endl;
    }

    pVectorsPtr_ = new edgeVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    edgeVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new edgeVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    edgeVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const areaVectorField& C = mesh_.areaCentres();
    const edgeScalarField& w = mesh_.weights();


    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nFaces(), symmTensor::zero);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        symmTensor wdd = (1.0/magSqr(d))*sqr(d);

        dd[own] += wdd;
        dd[nei] += wdd;
    }


    forAll(lsP.boundaryField(), patchi)
    {
        const faePatchScalarField& pw = w.boundaryField()[patchi];
        // Note: least squares in 1.4.1 and other OpenCFD versions
        // are wrong because of incorrect weighting.  HJ, 23/Oct/2008
//         const faePatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const faPatch& p = pw.patch();
        const unallocLabelList& edgeFaces = p.edgeFaces();

        // Build the d-vectors
        // HJ, reconsider deltas at the boundary, consistent with FVM
        // Current implementation is good for fixedValue boudaries, but may
        // cause problems with fixedGradient.  HJ, 4/Oct/2010
        vectorField pd = p.delta();

        if (p.coupled())
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                dd[edgeFaces[patchFacei]] +=
                    (1.0/magSqr(d))*sqr(d);
            }
        }
        else
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                dd[edgeFaces[patchFacei]] +=
                    (1.0/magSqr(d))*sqr(d);
            }
        }
    }


    // Invert the dd tensor
    symmTensorField invDd = inv(dd);


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        scalar magSfByMagSqrd = 1.0/magSqr(d);

        lsP[facei] = magSfByMagSqrd*(invDd[own] & d);
        lsN[facei] = -magSfByMagSqrd*(invDd[nei] & d);
    }

    forAll(lsP.boundaryField(), patchi)
    {
        faePatchVectorField& patchLsP = lsP.boundaryField()[patchi];

        const faePatchScalarField& pw = w.boundaryField()[patchi];
        // Note: least squares in 1.4.1 and other OpenCFD versions
        // are wrong becaus eof incorrect weighting.  HJ, 23/Oct/2008
//         const faePatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const faPatch& p = pw.patch();
        const unallocLabelList& edgeFaces = p.edgeFaces();

        // Build the d-vectors
        vectorField pd = p.delta();

        if (p.coupled())
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                patchLsP[patchFacei] =
                    (1.0/magSqr(d))
                   *(invDd[edgeFaces[patchFacei]] & d);
            }
        }
        else
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                patchLsP[patchFacei] =
                    (1.0/magSqr(d))
                   *(invDd[edgeFaces[patchFacei]] & d);
            }
        }
    }

    if (debug)
    {
        Info<< "leastSquaresFaVectors::makeLeastSquaresVectors() :"
            << "Finished constructing finite area least square gradient vectors"
            << endl;
    }
}


const Foam::edgeVectorField& Foam::leastSquaresFaVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::edgeVectorField& Foam::leastSquaresFaVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


bool Foam::leastSquaresFaVectors::movePoints() const
{
    if (debug)
    {
        InfoIn("bool leastSquaresFaVectors::movePoints() const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}

bool Foam::leastSquaresFaVectors::updateMesh(const mapPolyMesh&) const
{
    if (debug)
    {
        InfoIn("bool leastSquaresFaVectors::updateMesh() const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}

// ************************************************************************* //
