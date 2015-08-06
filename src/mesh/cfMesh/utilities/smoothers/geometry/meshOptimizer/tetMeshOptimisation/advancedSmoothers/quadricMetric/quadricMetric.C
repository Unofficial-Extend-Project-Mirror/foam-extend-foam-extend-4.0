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
#include "quadricMetric.H"
#include "partTetMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private member functions

scalar quadricMetric::evaluateMetric() const
{
    scalar val(0.0);

    forAll(normals_, nI)
        val += Foam::sqr(normals_[nI] & (p_ - centres_[nI]));

    return val;
}

void quadricMetric::evaluateGradients(vector& grad, tensor& gradGrad) const
{
    grad = vector::zero;
    gradGrad = tensor::zero;

    forAll(normals_, nI)
    {
        const scalar fx = (normals_[nI] & (p_ - centres_[nI]));

        grad += fx * normals_[nI];
        gradGrad += normals_[nI] * normals_[nI];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

quadricMetric::quadricMetric(partTetMeshSimplex& simplex)
:
    simplexSmoother(simplex),
    p_(simplex.pts()[simplex.tets()[0][3]]),
    normals_(),
    centres_()
{
    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const triangle<point, point> tri
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()]
        );

        const vector n = tri.normal();
        const scalar d = mag(n);

        if( d > VSMALL )
        {
            centres_.append(tri.centre());
            normals_.append(n/d);
        }
    }
}


quadricMetric::~quadricMetric()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation of knupp metric untangling
void quadricMetric::optimizeNodePosition(const scalar tolObsolete)
{
    if( !bb_.contains(p_) )
        p_ = 0.5 * (bb_.min() + bb_.max());

    const scalar tol = Foam::sqr(2.0 * SMALL) * magSqr(bb_.min() - bb_.max());

    label iterI, outerIter(0);

    vector gradF, disp;
    tensor gradGradF;

    scalar func, lastFunc;

    # ifdef DEBUGSmooth
    forAll(normals_, nI)
    {
        const scalar fx = normals_[nI] & (p_ - centres_[nI]);
        Info << "Tet " << nI << " has distance " << fx << " func "
            << Foam::sqr(mag(fx) - fx) << endl;
    }
    Info << "BoundBox size " << (bb_.max() - bb_.min()) << endl;
    Info << "Tolerance " << tol << endl;
    # endif

    bool finished;
    do
    {
        finished = true;

        lastFunc = evaluateMetric();

        iterI = 0;
        do
        {
            # ifdef DEBUGSmooth
            Info << "Iteration " << iterI << endl;
            Info << "Initial metric value " << lastFunc << endl;
            # endif

            //- store previous value
            const point pOrig = p_;

            //- evaluate gradients
            evaluateGradients(gradF, gradGradF);

            //- calculate displacement
            const scalar determinant = det(gradGradF);
            if( determinant > SMALL )
            {
                disp = (inv(gradGradF, determinant) & gradF);

                for(direction i=0;i<vector::nComponents;++i)
                {
                    const scalar& val = disp[i];
                    if( (val != val) || ((val - val) != (val - val)) )
                    {
                        disp = vector::zero;
                        break;
                    }
                }

                p_ -= disp;

                func = evaluateMetric();

                # ifdef DEBUGSmooth
                Info << "Second grad " << gradGradF << endl;
                Info << "inv(gradGradF, determinant) "
                    << inv(gradGradF, determinant) << endl;
                Info << "Gradient " << gradF << endl;
                Info << "Determinant " << determinant << endl;
                Info << "Displacement " << disp << endl;
                Info << "New metric value " << func << endl;
                # endif

                scalar relax(0.8);
                label nLoops(0);
                while( func > lastFunc )
                {
                    p_ = pOrig - relax * disp;
                    relax *= 0.5;
                    func = evaluateMetric();

                    if( func < lastFunc )
                        continue;

                    //- it seems that this direction is wrong
                    if( ++nLoops == 5 )
                    {
                        p_ = pOrig;
                        disp = vector::zero;
                        func = 0.0;
                    }
                }

                lastFunc = func;
            }
            else
            {
                disp = vector::zero;
            }
        } while( (magSqr(disp) > tol) && (++iterI < 10) );

        if( lastFunc < VSMALL )
            finished = false;
    } while( !finished && (++outerIter < 5) );

    # ifdef DEBUGSmooth
    Info << "Last value " << lastFunc << endl;
    Info << "Beta " << beta_ << endl;
    Info << "Metric with no beta " << evaluateMetricNoBeta() << endl;
    forAll(normals_, nI)
    {
        const scalar fx = normals_[nI] & (p_ - centres_[nI]);
        Info << "Tet " << nI << " has distance " << fx << " func "
            << Foam::sqr(mag(fx) - fx) << endl;
    }
    //::exit(1);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
