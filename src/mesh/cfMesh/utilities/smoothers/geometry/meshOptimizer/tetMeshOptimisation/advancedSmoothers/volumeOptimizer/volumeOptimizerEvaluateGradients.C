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
#include "volumeOptimizer.H"
#include "tetrahedron.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar volumeOptimizer::evaluateFunc() const
{
    const scalar K = evaluateStabilisationFactor();

    scalar func(0.0);

    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const tetrahedron<point, point> tet
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()],
            points_[pt.d()]
        );

        const scalar LSqrTri
        (
            magSqr(tet.d() - tet.a()) +
            magSqr(tet.d() - tet.b()) +
            magSqr(tet.d() - tet.c())
        );

        const scalar Vtri = tet.mag();
        const scalar stab = sqrt(sqr(Vtri) + K);
        const scalar Vstab = 0.5 * (Vtri + stab);

        func += LSqrTri / pow(Vstab, 2./3.);
    }

    return func;
}

scalar volumeOptimizer::evaluateStabilisationFactor() const
{
    scalar K = 0.0;

    scalar Vmin(VGREAT), LSqMax(0.0);

    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const tetrahedron<point, point> tet
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()],
            points_[pt.d()]
        );

        const scalar Vtri = tet.mag();

        Vmin = Foam::min(Vmin, Vtri);

        const scalar LSqrTri
        (
            magSqr(tet.d() - tet.a()) +
            magSqr(tet.d() - tet.b()) +
            magSqr(tet.d() - tet.c())
        );

        LSqMax = Foam::max(LSqMax, LSqrTri);
    }

    if( Vmin < SMALL * LSqMax )
        K = SMALL * LSqMax;

    return K;
}

void volumeOptimizer::evaluateGradientsExact
(
    vector& gradF,
    tensor& gradGradF
) const
{
    gradF = vector::zero;
    gradGradF = tensor::zero;

    const scalar K = evaluateStabilisationFactor();

    tensor gradGradLsq(tensor::zero);
    gradGradLsq.xx() = 6.0;
    gradGradLsq.yy() = 6.0;
    gradGradLsq.zz() = 6.0;

    const point& p = points_[pointI_];

    const scalar C = 2./3. * pow(0.5, 2./3.);
    const scalar C1 = C / 3.;

    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const tetrahedron<point, point> tet
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()],
            points_[pt.d()]
        );

        //- calculate the gradient of the volume
        const vector gradV
        (
            (1.0/6.0) *
            (
                (tet.b() - tet.a()) ^
                (tet.c() - tet.a())
            )
        );

        //- calculate the Frobenius norm
        const scalar LSqrTri
        (
            magSqr(tet.d() - tet.a()) +
            magSqr(tet.d() - tet.b()) +
            magSqr(tet.d() - tet.c())
        );

        //- calculate the volume of the tetrahedron
        const scalar Vtri = tet.mag();

        //- calculate the stabilisation factor for the volume
        const scalar stab = sqrt(sqr(Vtri) + K);

        //- evaluate the stabilised volume
        const scalar Vs = 0.5 * (Vtri + stab);

        if( Vs < VSMALL )
        {
            Info << "Tet " << tet << endl;
            Info << "gradV " << gradV << endl;
            Info << "Vtri " << Vtri << endl;
            IOstream::defaultPrecision(20);
            Info << "Vstab " << Vs << endl;

            FatalErrorIn
            (
                "void nodeDisplacementVolumeOptimizer()"
            ) << "I cannot continue " << exit(FatalError);
        }

        //- calculate the gradient of the stabilisation volume
        const vector gradStab = Vtri * gradV / stab;

        //- calculate the gradient of the Frobenius norm
        const vector gradLsq = 2. * (3. * p - tet.a() - tet.b() - tet.c());

        //- calculate the gradient of the stabilised volume
        const vector gradVs = 0.5 * (gradV + gradStab);

        //- calculate the gradient of the functional
        const scalar Vs13 = pow(2. * Vs, 1./3.);
        const scalar Vstab = pow(Vs, 2./3.);
        const scalar sqrVstab = sqr(Vstab);
        const vector gradVstab = C * (2. * gradVs) / Vs13;

        gradF += gradLsq / Vstab - LSqrTri * gradVstab / sqrVstab;

        //- calculate the second gradient of the stabilisation volume
        const tensor gradGradStab =
            (gradV * gradV) / stab -
            sqr(Vtri) * (gradV * gradV) / pow(stab, 3);

        //- calculate the second gradient of the stabilised volume

        const tensor gradGradVstab =
            C * (gradGradStab / Vs13) -
            C1 * 4. * (gradVs * gradVs) / pow(Vs13, 4);

        //- calculate the second gradient
        gradGradF +=
            gradGradLsq / Vstab -
            twoSymm(gradLsq * (gradVstab / sqrVstab)) -
            LSqrTri * gradGradVstab / sqrVstab +
            2.0 * LSqrTri * (gradVstab * gradVstab) / (sqrVstab * Vstab);
    }
}

scalar volumeOptimizer::optimiseDivideAndConquer(const scalar tol)
{
    point& pOpt = points_[pointI_];

    pOpt = 0.5 * (bb_.max() + bb_.min());
    point currCentre = pOpt;
    scalar dx = (bb_.max().x() - bb_.min().x()) / 2.0;
    scalar dy = (bb_.max().y() - bb_.min().y()) / 2.0;
    scalar dz = (bb_.max().z() - bb_.min().z()) / 2.0;

    FixedList<vector, 8> dirVecs;
    dirVecs[0] = vector(-1.0, -1.0, -1.0);
    dirVecs[1] = vector(1.0, -1.0, -1.0);
    dirVecs[2] = vector(-1.0, 1.0, -1.0);
    dirVecs[3] = vector(1.0, 1.0, -1.0);
    dirVecs[4] = vector(-1.0, -1.0, 1.0);
    dirVecs[5] = vector(1.0, -1.0, 1.0);
    dirVecs[6] = vector(-1.0, 1.0, 1.0);
    dirVecs[7] = vector(1.0, 1.0, 1.0);

    label iter(0);

    //- find the value of the functional in the centre of the bnd box
    scalar funcBefore, funcAfter(evaluateFunc());

    do
    {
        funcBefore = funcAfter;

        funcAfter = VGREAT;
        point minCentre(vector::zero);

        forAll(dirVecs, i)
        {
            pOpt.x() = currCentre.x() + 0.5 * dirVecs[i].x() * dx;
            pOpt.y() = currCentre.y() + 0.5 * dirVecs[i].y() * dy;
            pOpt.z() = currCentre.z() + 0.5 * dirVecs[i].z() * dz;

            const scalar func = evaluateFunc();

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
        dz *= 0.5;

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

scalar volumeOptimizer::optimiseSteepestDescent(const scalar tol)
{
    label iter(0);

    point& p = points_[pointI_];

    # ifdef DEBUGSmooth
    Info << nl << "Smoothing point " << pointI_
         << " with coordinates " << p << endl;
    scalar Vmina(VGREAT);
    forAll(tets_, tetI)
    Vmina = Foam::min(Vmina, tets_[tetI].mag(points_));
    Info << "Vmin before " << Vmina << endl;
    # endif

    vector gradF;
    vector disp(vector::zero);
    tensor gradGradF;
    point pOrig;

    scalar funcBefore, funcAfter(evaluateFunc());

    bool finished;
    do
    {
        finished = false;
        pOrig = p;
        funcBefore = funcAfter;

        evaluateGradientsExact(gradF, gradGradF);

        const scalar determinant = Foam::det(gradGradF);
        if( determinant > SMALL )
        {
            disp = (inv(gradGradF, determinant) & gradF);

            p -= disp;

            funcAfter = evaluateFunc();

            # ifdef DEBUGSmooth
            Info << nl << "gradF " << gradF << endl;
            Info << "gradGradF " << gradGradF << endl;
            Info << "det(gradGradF) " << determinant << endl;
            Info << "disp " << disp << endl;
            Info << "Func before " << funcBefore << endl;
            Info << "Func after " << funcAfter << endl;
            # endif

            scalar relax(0.8);
            label nLoops(0);

            while( funcAfter > funcBefore )
            {
                p = pOrig - relax * disp;
                relax *= 0.5;
                funcAfter = evaluateFunc();

                if( funcAfter < funcBefore )
                    continue;

                if( ++nLoops == 5 )
                {
                    //- it seems that this direction is wrong, stop the loop
                    p = pOrig;
                    disp = vector::zero;
                    finished = true;
                    funcAfter = funcBefore;
                }
            }

            if( mag(funcBefore - funcAfter) / funcBefore < tol )
                finished = true;
        }
        else
        {
            //- move in random direction
            //- this is usually needed to move the point off the zero volume
            disp = vector::zero;
            forAll(tets_, tetI)
            {
                const partTet& tet = tets_[tetI];
                const scalar Vtri = tet.mag(points_);

                if( Vtri < SMALL )
                {
                    triangle<point, point> tri
                    (
                        points_[tet.a()],
                        points_[tet.b()],
                        points_[tet.c()]
                    );

                    vector n = tri.normal();
                    const scalar d = mag(n);

                    if( d > VSMALL )
                        disp += 0.01 * (n / d);
                }
            }

            p += disp;
            funcAfter = evaluateFunc();
        }
    } while( (++iter < 100) && !finished );

    # ifdef DEBUGSmooth
    scalar Vmin(VGREAT);
    forAll(tets_, tetI)
        Vmin = Foam::min(Vmin, tets_[tetI].mag(points_));

    Info << nl << "New coordinates for point "
        << pointI_ << " are " << p << endl;
    Info << "Num iterations " << iter << " gradient " << gradF << endl;
    Info << "Vmin " << Vmin << endl;
    # endif

    return funcAfter;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
