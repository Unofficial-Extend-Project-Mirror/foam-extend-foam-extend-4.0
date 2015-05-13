/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"
#include "surfaceOptimizer.H"
#include "matrix2D.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar surfaceOptimizer::evaluateStabilisationFactor() const
{
    //- find the minimum area
    scalar Amin(VGREAT), LsqMax(0.0);
    forAll(trias_, triI)
    {
        const point& p0 = pts_[trias_[triI][0]];
        const point& p1 = pts_[trias_[triI][1]];
        const point& p2 = pts_[trias_[triI][2]];

        const scalar Atri =
            0.5 *
            (
                (p1.x() - p0.x()) * (p2.y() - p0.y()) -
                (p2.x() - p0.x()) * (p1.y() - p0.y())
            );

        # ifdef DEBUGSmooth
        Info << "Triangle " << triI << " area " << Atri << endl;
        # endif

        const scalar LSqrTri
        (
            magSqr(p0 - p1) +
            magSqr(p2 - p0)
        );

        Amin = Foam::min(Amin, Atri);
        LsqMax = Foam::max(LsqMax, LSqrTri);
    }

    # ifdef DEBUGSmooth
    Info << "Amin " << Amin << endl;
    Info << "LsqMax " << LsqMax << endl;
    # endif

    //- K is greater than zero in case the stabilisation is needed
    scalar K = 0.0;
    if( Amin < SMALL * LsqMax )
    {
        K = SMALL * LsqMax;
    }

    return K;
}

scalar surfaceOptimizer::evaluateFunc(const scalar& K) const
{
    scalar func(0.0);

    forAll(trias_, triI)
    {
        const point& p0 = pts_[trias_[triI][0]];
        const point& p1 = pts_[trias_[triI][1]];
        const point& p2 = pts_[trias_[triI][2]];

        const scalar Atri =
            0.5 *
            (
                (p1.x() - p0.x()) * (p2.y() - p0.y()) -
                (p2.x() - p0.x()) * (p1.y() - p0.y())
            );

        const scalar stab = sqrt(sqr(Atri) + K);

        # ifdef DEBUGSmooth
        Info << "Triangle " << triI << " area " << Atri << endl;
        # endif

        const scalar LSqrTri
        (
            magSqr(p0 - p1) +
            magSqr(p2 - p0)
        );

        func += LSqrTri / (0.5 * (Atri + stab));
    }

    return func;
}

//- evaluate gradients needed for optimisation
void surfaceOptimizer::evaluateGradients
(
    const scalar& K,
    vector& gradF,
    tensor& gradGradF
) const
{
    gradF = vector::zero;
    gradGradF = tensor::zero;

    tensor gradGradLt(tensor::zero);
    gradGradLt.xx() = 4.0;
    gradGradLt.yy() = 4.0;

    forAll(trias_, triI)
    {
        const point& p0 = pts_[trias_[triI][0]];
        const point& p1 = pts_[trias_[triI][1]];
        const point& p2 = pts_[trias_[triI][2]];

        if( magSqr(p1 - p2) < VSMALL ) continue;

        const scalar LSqrTri
        (
            magSqr(p0 - p1) +
            magSqr(p2 - p0)
        );

        const scalar Atri =
            0.5 *
            (
                (p1.x() - p0.x()) * (p2.y() - p0.y()) -
                (p2.x() - p0.x()) * (p1.y() - p0.y())
            );

        const scalar stab = sqrt(sqr(Atri) + K);

        const scalar Astab = 0.5 * (Atri + stab);

        const vector gradAtri
        (
            0.5 * (p1.y() - p2.y()),
            0.5 * (p2.x() - p1.x()),
            0.0
        );

        const vector gradAstab = 0.5 * (gradAtri + Atri * gradAtri / stab);

        const tensor gradGradAstab =
            0.5 *
            (
                (gradAtri * gradAtri) / stab -
                sqr(Atri) * (gradAtri * gradAtri) / pow(stab, 3)
            );

        const vector gradLt(4.0 * p0 - 2.0 * p1 - 2.0 * p2);

        //- calculate the gradient
        const scalar sqrAstab = sqr(Astab);
        gradF += gradLt / Astab - (LSqrTri * gradAstab) / sqrAstab;

        //- calculate the second gradient
        gradGradF +=
            gradGradLt / Astab -
            twoSymm(gradLt * gradAstab) / sqrAstab -
            gradGradAstab * LSqrTri / sqrAstab +
            2.0 * LSqrTri * (gradAstab * gradAstab) / (sqrAstab * Astab);
    }

    //- stabilise diagonal terms
    if( mag(gradGradF.xx()) < VSMALL ) gradGradF.xx() = VSMALL;
    if( mag(gradGradF.yy()) < VSMALL ) gradGradF.yy() = VSMALL;
}

scalar surfaceOptimizer::optimiseDivideAndConquer(const scalar tol)
{
    point& pOpt = pts_[trias_[0][0]];

    pOpt = 0.5 * (pMax_ + pMin_);
    point currCentre = pOpt;
    scalar dx = (pMax_.x() - pMin_.x()) / 2.0;
    scalar dy = (pMax_.y() - pMin_.y()) / 2.0;

    FixedList<vector, 4> dirVecs;
    dirVecs[0] = vector(-1.0, -1.0, 0.0);
    dirVecs[1] = vector(1.0, -1.0, 0.0);
    dirVecs[2] = vector(-1.0, 1.0, 0.0);
    dirVecs[3] = vector(1.0, 1.0, 0.0);

    label iter(0);

    //- find the value of the functional in the centre of the bnd box
    scalar K = evaluateStabilisationFactor();
    scalar funcBefore, funcAfter(evaluateFunc(K));

    do
    {
        funcBefore = funcAfter;

        funcAfter = VGREAT;
        point minCentre(vector::zero);

        forAll(dirVecs, i)
        {
            pOpt.x() = currCentre.x() + 0.5 * dirVecs[i].x() * dx;
            pOpt.y() = currCentre.y() + 0.5 * dirVecs[i].y() * dy;

            K = evaluateStabilisationFactor();
            const scalar func = evaluateFunc(K);

            if( func < funcAfter )
            {
                minCentre = pOpt;
                funcAfter = func;
            }
        }

        //- set the centre with the minimum value
        //- as the centre for future search
        currCentre = minCentre;
        pOpt = minCentre;

        //- halve the search range
        dx *= 0.5;
        dy *= 0.5;

        //- calculate the tolerence
        const scalar t = mag(funcAfter - funcBefore) / funcAfter;

        # ifdef DEBUGSmooth
        Info << "Point position " << pOpt << endl;
        Info << "Func before " << funcBefore << endl;
        Info << "Func after " << funcAfter << endl;
        Info << "Normalised difference " << t << endl;
        # endif

        if( t < tol )
            break;
    } while( ++iter < 100 );

    return funcAfter;
}

scalar surfaceOptimizer::optimiseSteepestDescent(const scalar tol)
{
    point& pOpt = pts_[trias_[0][0]];

    //- find the bounding box
    const scalar avgEdge = Foam::mag(pMax_ - pMin_);

    //- find the minimum value on the 5 x 5 raster
    scalar K = evaluateStabilisationFactor();
    scalar funcBefore, funcAfter(evaluateFunc(K));

    //- start with steepest descent optimisation
    vector gradF;
    tensor gradGradF;
    vector disp;
    disp.z() = 0.0;

    direction nIterations(0);
    do
    {
        funcBefore = funcAfter;

        evaluateGradients(K, gradF, gradGradF);

        //- store data into a matrix
        matrix2D mat;
        mat[0][0] = gradGradF.xx();
        mat[0][1] = gradGradF.xy();
        mat[1][0] = gradGradF.yx();
        mat[1][1] = gradGradF.yy();
        FixedList<scalar, 2> source;
        source[0] = gradF.x();
        source[1] = gradF.y();

        //- calculate the determinant
        const scalar det = mat.determinant();

        if( mag(det) < VSMALL )
        {
            disp = vector::zero;
        }
        else
        {
            disp.x() = mat.solveFirst(source);
            disp.y() = mat.solveSecond(source);

            if( mag(disp) > 0.2 * avgEdge )
            {
                vector dir = disp / mag(disp);

                disp = dir * 0.2 * avgEdge;
            }
        }

        # ifdef DEBUGSmooth
        Info << "Second gradient " << gradGradF << endl;
        Info << "Gradient " << gradF << endl;
        Info << "Displacement " << disp << endl;
        Info << "K = " << K << endl;
        # endif

        pOpt -= disp;

        K = evaluateStabilisationFactor();
        funcAfter = evaluateFunc(K);

        if( mag(funcAfter - funcBefore) / funcBefore < tol )
            break;

        #ifdef DEBUGSmooth
        Info << "New coordinates " << pOpt << endl;
        # endif

    } while( ++nIterations < 100 );

    return funcAfter;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

surfaceOptimizer::surfaceOptimizer
(
    DynList<point>& pts,
    const DynList<triFace>& trias
)
:
    pts_(pts),
    trias_(trias),
    pMin_(),
    pMax_()
{
    pMin_ = pts_[trias_[0][1]];
    pMax_ = pMin_;

    forAll(trias_, triI)
    {
        const triFace& tf = trias_[triI];

        for(label i=1;i<3;++i)
        {
            pMin_ = Foam::min(pMin_, pts_[tf[i]]);
            pMax_ = Foam::max(pMax_, pts_[tf[i]]);
        }
    }
}

surfaceOptimizer::~surfaceOptimizer()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point surfaceOptimizer::optimizePoint(const scalar tol)
{
    const scalar scale = mag(pMax_ - pMin_);
    forAll(pts_, i)
        pts_[i] /= scale;
    pMin_ /= scale;
    pMax_ /= scale;

    point& pOpt = pts_[trias_[0][0]];

    const scalar funcDivide = optimiseDivideAndConquer(tol);
    const point newPoint = pOpt;

    const scalar funcSteepest = optimiseSteepestDescent(tol);

    if( funcSteepest > funcDivide )
        pOpt = newPoint;

    forAll(pts_, i)
        pts_[i] *= scale;
    pMin_ *= scale;
    pMax_ *= scale;

    return pOpt;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
