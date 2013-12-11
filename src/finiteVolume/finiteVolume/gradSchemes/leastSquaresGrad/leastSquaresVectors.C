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
    const surfaceScalarField& w = mesh().weights();
//     const surfaceScalarField& magSf = mesh().magSf();


    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh().nCells(), symmTensor::zero);

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
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        // Note: least squares in 1.4.1 and other OpenCFD versions
        // are wrong because of incorrect weighting.  HJ, 23/Oct/2008
//         const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const unallocLabelList& faceCells = p.patch().faceCells();

        // Build the d-vectors

        // Original version: closest distance to boundary
        vectorField pd =
            mesh().Sf().boundaryField()[patchi]
           /(
               mesh().magSf().boundaryField()[patchi]
              *mesh().deltaCoeffs().boundaryField()[patchi]
           );

        if (!mesh().orthogonal())
        {
            pd -= mesh().correctionVectors().boundaryField()[patchi]
                /mesh().deltaCoeffs().boundaryField()[patchi];
        }

        // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
        // Experimental: review fixed gradient condition.  HJ, 30/Sep/2010
//         vectorField pd = p.delta();

        if (p.coupled())
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                dd[faceCells[patchFacei]] += (1.0/magSqr(d))*sqr(d);
            }
        }
        else
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                dd[faceCells[patchFacei]] += (1.0/magSqr(d))*sqr(d);
            }
        }
    }


    // Invert the dd tensor
//     symmTensorField invDd = inv(dd);
    // Fix: householder inverse.  HJ, 3/Nov/2009
    symmTensorField invDd = hinv(dd);


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
        fvsPatchVectorField& patchLsP = lsP.boundaryField()[patchi];

        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        // Note: least squares in 1.4.1 and other OpenCFD versions
        // are wrong because of incorrect weighting.  HJ, 23/Oct/2008
//         const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const unallocLabelList& faceCells = p.faceCells();

        // Build the d-vectors
        // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
        vectorField pd = p.delta();

        if (p.coupled())
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                patchLsP[patchFacei] =
                    (1.0/magSqr(d))
                   *(invDd[faceCells[patchFacei]] & d);
            }
        }
        else
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                patchLsP[patchFacei] =
                    (1.0/magSqr(d))
                   *(invDd[faceCells[patchFacei]] & d);
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
        InfoIn("bool leastSquaresVectors::updateMesh(const mapPolyMesh&) const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}

// ************************************************************************* //
