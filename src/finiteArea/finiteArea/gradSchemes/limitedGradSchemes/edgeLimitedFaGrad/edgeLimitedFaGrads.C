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

#include "edgeLimitedFaGrad.H"
#include "gaussFaGrad.H"
#include "faMesh.H"
#include "areaMesh.H"
#include "edgeMesh.H"
#include "areaFields.H"
#include "fixedValueFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFaGradScheme(edgeLimitedGrad)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void edgeLimitedGrad<Type>::limitEdge
(
    scalar& limiter,
    const scalar maxDelta,
    const scalar minDelta,
    const scalar extrapolate
) const
{
    if (extrapolate > maxDelta + VSMALL)
    {
        limiter = min(limiter, maxDelta/extrapolate);
    }
    else if (extrapolate < minDelta - VSMALL)
    {
        limiter = min(limiter, minDelta/extrapolate);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<areaVectorField> edgeLimitedGrad<scalar>::grad
(
    const areaScalarField& vsf
) const
{
    const faMesh& mesh = vsf.mesh();

    tmp<areaVectorField> tGrad = basicGradScheme_().grad(vsf);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    areaVectorField& g = tGrad();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const areaVectorField& C = mesh.areaCentres();
    const edgeVectorField& Cf = mesh.edgeCentres();

    // create limiter
    scalarField limiter(vsf.internalField().size(), 1.0);

    scalar rk = (1.0/k_ - 1.0);

    forAll(owner, edgei)
    {
        label own = owner[edgei];
        label nei = neighbour[edgei];

        scalar vsfOwn = vsf[own];
        scalar vsfNei = vsf[nei];

        scalar maxEdge = max(vsfOwn, vsfNei);
        scalar minEdge = min(vsfOwn, vsfNei);
        scalar maxMinEdge = rk*(maxEdge - minEdge);
        maxEdge += maxMinEdge;
        minEdge -= maxMinEdge;

        // owner side
        limitEdge
        (
            limiter[own],
            maxEdge - vsfOwn, minEdge - vsfOwn,
            (Cf[edgei] - C[own]) & g[own]
        );

        // neighbour side
        limitEdge
        (
            limiter[nei],
            maxEdge - vsfNei, minEdge - vsfNei,
            (Cf[edgei] - C[nei]) & g[nei]
        );
    }

    const areaScalarField::GeometricBoundaryField& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const faPatchScalarField& psf = bsf[patchi];

        const unallocLabelList& pOwner = mesh.boundary()[patchi].edgeFaces();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        if (psf.coupled())
        {
            scalarField psfNei = psf.patchNeighbourField();

            forAll(pOwner, pEdgei)
            {
                label own = pOwner[pEdgei];

                scalar vsfOwn = vsf[own];
                scalar vsfNei = psfNei[pEdgei];

                scalar maxEdge = max(vsfOwn, vsfNei);
                scalar minEdge = min(vsfOwn, vsfNei);
                scalar maxMinEdge = rk*(maxEdge - minEdge);
                maxEdge += maxMinEdge;
                minEdge -= maxMinEdge;

                limitEdge
                (
                    limiter[own],
                    maxEdge - vsfOwn, minEdge - vsfOwn,
                    (pCf[pEdgei] - C[own]) & g[own]
                );
            }
        }
        else if (psf.fixesValue())
        {
            forAll(pOwner, pEdgei)
            {
                label own = pOwner[pEdgei];

                scalar vsfOwn = vsf[own];
                scalar vsfNei = psf[pEdgei];

                scalar maxEdge = max(vsfOwn, vsfNei);
                scalar minEdge = min(vsfOwn, vsfNei);
                scalar maxMinEdge = rk*(maxEdge - minEdge);
                maxEdge += maxMinEdge;
                minEdge -= maxMinEdge;

                limitEdge
                (
                    limiter[own],
                    maxEdge - vsfOwn, minEdge - vsfOwn,
                    (pCf[pEdgei] - C[own]) & g[own]
                );
            }
        }
    }

    if (fa::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.internalField() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


template<>
tmp<areaTensorField> edgeLimitedGrad<vector>::grad
(
    const areaVectorField& vvf
) const
{
    const faMesh& mesh = vvf.mesh();

    tmp<areaTensorField> tGrad = basicGradScheme_().grad(vvf);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    areaTensorField& g = tGrad();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const areaVectorField& C = mesh.areaCentres();
    const edgeVectorField& Cf = mesh.edgeCentres();

    // create limiter
    scalarField limiter(vvf.internalField().size(), 1.0);

    scalar rk = (1.0/k_ - 1.0);

    forAll(owner, edgei)
    {
        label own = owner[edgei];
        label nei = neighbour[edgei];

        vector vvfOwn = vvf[own];
        vector vvfNei = vvf[nei];

        // owner side
        vector gradf = (Cf[edgei] - C[own]) & g[own];

        scalar vsfOwn = gradf & vvfOwn;
        scalar vsfNei = gradf & vvfNei;

        scalar maxEdge = max(vsfOwn, vsfNei);
        scalar minEdge = min(vsfOwn, vsfNei);
        scalar maxMinEdge = rk*(maxEdge - minEdge);
        maxEdge += maxMinEdge;
        minEdge -= maxMinEdge;

        limitEdge
        (
            limiter[own],
            maxEdge - vsfOwn, minEdge - vsfOwn,
            magSqr(gradf)
        );


        // neighbour side
        gradf = (Cf[edgei] - C[nei]) & g[nei];

        vsfOwn = gradf & vvfOwn;
        vsfNei = gradf & vvfNei;

        maxEdge = max(vsfOwn, vsfNei);
        minEdge = min(vsfOwn, vsfNei);

        limitEdge
        (
            limiter[nei],
            maxEdge - vsfNei, minEdge - vsfNei,
            magSqr(gradf)
        );
    }


    const areaVectorField::GeometricBoundaryField& bvf = vvf.boundaryField();

    forAll(bvf, patchi)
    {
        const faPatchVectorField& psf = bvf[patchi];

        const unallocLabelList& pOwner = mesh.boundary()[patchi].edgeFaces();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        if (psf.coupled())
        {
            vectorField psfNei = psf.patchNeighbourField();

            forAll(pOwner, pEdgei)
            {
                label own = pOwner[pEdgei];

                vector vvfOwn = vvf[own];
                vector vvfNei = psfNei[pEdgei];

                vector gradf = (pCf[pEdgei] - C[own]) & g[own];

                scalar vsfOwn = gradf & vvfOwn;
                scalar vsfNei = gradf & vvfNei;

                scalar maxEdge = max(vsfOwn, vsfNei);
                scalar minEdge = min(vsfOwn, vsfNei);
                scalar maxMinEdge = rk*(maxEdge - minEdge);
                maxEdge += maxMinEdge;
                minEdge -= maxMinEdge;

                limitEdge
                (
                    limiter[own],
                    maxEdge - vsfOwn, minEdge - vsfOwn,
                    magSqr(gradf)
                );
            }
        }
        else if (psf.fixesValue())
        {
            forAll(pOwner, pEdgei)
            {
                label own = pOwner[pEdgei];

                vector vvfOwn = vvf[own];
                vector vvfNei = psf[pEdgei];

                vector gradf = (pCf[pEdgei] - C[own]) & g[own];

                scalar vsfOwn = gradf & vvfOwn;
                scalar vsfNei = gradf & vvfNei;

                scalar maxEdge = max(vsfOwn, vsfNei);
                scalar minEdge = min(vsfOwn, vsfNei);
                scalar maxMinEdge = rk*(maxEdge - minEdge);
                maxEdge += maxMinEdge;
                minEdge -= maxMinEdge;

                limitEdge
                (
                    limiter[own],
                    maxEdge - vsfOwn, minEdge - vsfOwn,
                    magSqr(gradf)
                );
            }
        }
    }

    if (fa::debug)
    {
        Info<< "gradient limiter for: " << vvf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.internalField() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vvf, g);

    return tGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
