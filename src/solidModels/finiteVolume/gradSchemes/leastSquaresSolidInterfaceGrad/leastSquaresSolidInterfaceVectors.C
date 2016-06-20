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

#include "leastSquaresSolidInterfaceVectors.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "mapPolyMesh.H"
#include "solidInterface.H"
#include "IOReferencer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresSolidInterfaceVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresSolidInterfaceVectors::leastSquaresSolidInterfaceVectors
(const fvMesh& mesh)
:
    MeshObject<fvMesh, leastSquaresSolidInterfaceVectors>(mesh),
    pVectorsPtr_(NULL),
    nVectorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresSolidInterfaceVectors::~leastSquaresSolidInterfaceVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresSolidInterfaceVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "leastSquaresSolidInterfaceVectors::makeLeastSquaresVectors() :"
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
    //const surfaceScalarField& magSf = mesh().magSf();
    const surfaceVectorField& Cf = mesh().Cf();

    // interface fields
    const solidInterface& solInt =
      mesh().objectRegistry::lookupObject<IOReferencer<solidInterface> >
        ("solidInterface")();
    const labelList& interfaceFacesMap = solInt.indicatorFieldMap();
    const labelListList& interfaceProcPatchFacesMap =
        solInt.processorPatchFacesMap();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh().nCells(), symmTensor::zero);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        if (interfaceFacesMap[facei] > -SMALL)
        {
            // for interface faces, we use the face centre value
            // instead of the neighbour cell centre
            vector ownD = Cf[facei] - C[own];
            vector neiD = C[nei] - Cf[facei];
            symmTensor ownWdd = (1.0/magSqr(ownD))*sqr(ownD);
            symmTensor neiWdd = (1.0/magSqr(neiD))*sqr(neiD);

            dd[own] += ownWdd;
            dd[nei] += neiWdd;
        }
        else
        {
            // standard method
            vector d = C[nei] - C[own];
            symmTensor wdd = (1.0/magSqr(d))*sqr(d);

            dd[own] += wdd;
            dd[nei] += wdd;
        }
    }


    forAll(lsP.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        // Note: least squares in 1.4.1 and other OpenCFD versions
        // are wrong because of incorrect weighting.  HJ, 23/Oct/2008
        //         const fvsPatchScalarField& pMagSf =
        //magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const unallocLabelList& faceCells = p.patch().faceCells();

        // Build the d-vectors

        // Original version: closest distance to boundary
        // vectorField pd =
//             mesh().Sf().boundaryField()[patchi]
//            /(
//                mesh().magSf().boundaryField()[patchi]
//               *mesh().deltaCoeffs().boundaryField()[patchi]
//            );

//         if (!mesh().orthogonal())
//         {
//             pd -= mesh().correctionVectors().boundaryField()[patchi]
//                 /mesh().deltaCoeffs().boundaryField()[patchi];
//          }

        // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
        // Experimental: review fixed gradient condition.  HJ, 30/Sep/2010
        // philipc - full delta vector imperative for solid mechanics
        vectorField pd = p.delta();

        if (p.coupled())
        {
            forAll(pd, patchFacei)
            {
                // philipc: special treatment of solid inter faces
                if (interfaceProcPatchFacesMap[patchi][patchFacei] > -SMALL)
                {
                    // use face centre instead of neighbour
                    const vector d =
                        Cf.boundaryField()[patchi][patchFacei]
                        - C[faceCells[patchFacei]];

                    dd[faceCells[patchFacei]] += (1.0/magSqr(d))*sqr(d);
                }
                else
                {
                    // standard method
                    const vector& d = pd[patchFacei];

                    dd[faceCells[patchFacei]] += (1.0/magSqr(d))*sqr(d);
                }
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

        if (interfaceFacesMap[facei] > -SMALL)
        {
            // for interface faces, we use the face centre value
            // instead of the neighbour cell centre
            vector ownD = Cf[facei] - C[own];
            vector neiD = C[nei] - Cf[facei];
            scalar ownMagSfByMagSqrd = 1.0/magSqr(ownD);
            scalar neiMagSfByMagSqrd = 1.0/magSqr(neiD);

            lsP[facei] = ownMagSfByMagSqrd*(invDd[own] & ownD);
            lsN[facei] = -neiMagSfByMagSqrd*(invDd[nei] & neiD);
        }
        else
        {
            // standard method
            vector d = C[nei] - C[own];
            scalar magSfByMagSqrd = 1.0/magSqr(d);

            lsP[facei] = magSfByMagSqrd*(invDd[own] & d);
            lsN[facei] = -magSfByMagSqrd*(invDd[nei] & d);
        }
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
                // philipc: special treatment of solid inter faces
                if (interfaceProcPatchFacesMap[patchi][patchFacei] > -SMALL)
                {
                    // use face centre instead of neighbour
                    const vector d =
                        Cf.boundaryField()[patchi][patchFacei]
                        - C[faceCells[patchFacei]];

                    patchLsP[patchFacei] =
                        (1.0/magSqr(d))
                        *(invDd[faceCells[patchFacei]] & d);
                }
                else
                {
                    // standard method
                    const vector& d = pd[patchFacei];

                    patchLsP[patchFacei] =
                        (1.0/magSqr(d))
                        *(invDd[faceCells[patchFacei]] & d);
                }
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
        Info<< "leastSquaresSolidInterfaceVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::surfaceVectorField&
Foam::leastSquaresSolidInterfaceVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField&
Foam::leastSquaresSolidInterfaceVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


bool Foam::leastSquaresSolidInterfaceVectors::movePoints() const
{
    if (debug)
    {
        InfoIn("bool leastSquaresSolidInterfaceVectors::movePoints() const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}

bool Foam::leastSquaresSolidInterfaceVectors::updateMesh
(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn("bool leastSquaresSolidInterfaceVectors::updateMesh"
               "(const mapPolyMesh&) const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}

// ************************************************************************* //
