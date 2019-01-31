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

#include "quadricFitting.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * * * //

void quadricFitting::calculateNormalVector()
{
    symmTensor mat(symmTensor::zero);

    forAll(otherPoints_, i)
    {
        const vector d = (otherPoints_[i] - origin_);

        mat += symm(d * d);
    }

    mat /= otherPoints_.size();

    //- find the eigenvalues of the tensor
    const vector ev = help::eigenValues(mat);

    //- estimate the normal as the eigenvector associated
    //- to the smallest eigenvalue
    # ifdef OpenCFDSpecific
    normal_ = eigenVectors(mat, ev).x();
    # else
    normal_ = help::eigenVector(mat, ev[0]);
    # endif
}

void quadricFitting::calculateCoordinateSystem()
{
    if( mag(normal_) < VSMALL )
    {
        vecX_ = vecY_ = vector::zero;
        forAll(transformedPoints_, i)
            transformedPoints_[i] = vector::zero;

        return;
    }

    const plane pl(origin_, normal_);

    bool valid(false);

    label pI(0);
    while( !valid )
    {
        const point vp = pl.nearestPoint(otherPoints_[pI]);
        vecX_ = vp - origin_;

        if( magSqr(vecX_) < VSMALL )
        {
            ++pI;
            continue;
        }

        vecX_ /= mag(vecX_);
        vecY_ = normal_ ^ vecX_;
        valid = true;
    }

    if( !valid )
    {
        FatalErrorIn
        (
            "void quadricFitting::calculateCoordinateSystem()"
        ) << "Cannot find a valid coordinate system!" << endl;

        normal_ = vecX_ = vecY_ = vector::zero;
    }

    //- transform other points to this new coordinate system
    transformedPoints_.setSize(otherPoints_.size());

    forAll(otherPoints_, i)
    {
        const vector delta = (otherPoints_[i] - origin_);
        point& p = transformedPoints_[i];

        p.x() = delta & vecX_;
        p.y() = delta & vecY_;
        p.z() = delta & normal_;
    }
}

void quadricFitting::calculateQuadricCoeffs()
{
    if( mag(normal_) < VSMALL )
    {
        coefficients_.setSize(5);
        coefficients_ = 0.0;
        return;
    }

    simpleMatrix<scalar> mat(5, 0.0, 0.0);
    DynList<scalar> helper;
    helper.setSize(5);

    forAll(transformedPoints_, i)
    {
        const point& p = transformedPoints_[i];

        helper[0] = sqr(p.x());
        helper[1] = sqr(p.y());
        helper[2] = p.x() * p.y();
        helper[3] = p.x();
        helper[4] = p.y();

        for(label rowI=0;rowI<5;++rowI)
        {
            for(label colI=rowI;colI<5;++colI)
                mat[rowI][colI] += helper[rowI] * helper[colI];

            mat.source()[rowI] += helper[rowI] * p.z();
        }
    }

    for(label rowI=1;rowI<5;++rowI)
    {
        for(label colI=0;colI<rowI;++colI)
            mat[rowI][colI] = mat[colI][rowI];

        if( mag(mat[rowI][rowI]) < SMALL )
            mat[rowI][rowI] = SMALL;
    }

    coefficients_ = mat.LUsolve();
}

void quadricFitting::calculateBestFit()
{
    label iteration(0);

    while( iteration++ < 10 )
    {
        calculateCoordinateSystem();

        calculateQuadricCoeffs();

        if( (mag(coefficients_[3]) > SMALL) || (mag(coefficients_[4]) > SMALL) )
        {
            //- correct the normal
            const scalar d =
                Foam::sqrt(1 + sqr(coefficients_[3]) + sqr(coefficients_[4]));

            const vector newNormal
            (
                normal_ / d -
                coefficients_[3] * vecX_ / d -
                coefficients_[4] * vecY_ / d
            );

            normal_ = newNormal;
        }
        else
        {
            break;
        }
    }

    # ifdef DEBUGQuadric
    Info << "Other points " << otherPoints_ << endl;
    Info << "Transformed points " << transformedPoints_ << endl;
    Info << "normal " << normal_ << endl;
    Info << "Coefficients " << coefficients_ << endl;
    # endif
}

// * * * * * * * * * * * * * * * *  Constructors * * * * * * * * * * * * * * //

template<class ListType>
inline quadricFitting::quadricFitting
(
    const point& origin,
    const vector normal,
    const ListType& otherPoints
)
:
    origin_(origin),
    normal_(normal),
    vecX_(vector::zero),
    vecY_(vector::zero),
    otherPoints_(),
    transformedPoints_(),
    coefficients_()
{
    otherPoints_.setSize(otherPoints.size());
    forAll(otherPoints, i)
        otherPoints_[i] = otherPoints[i];

    if( magSqr(normal_) < VSMALL )
    {
        FatalErrorIn
        (
            "template<class ListType>\n"
            "inline quadricFitting::quadricFitting"
            "(const point& origin, const vector normal,"
            " const ListType& otherPoints)"
        ) << "Cannot construct the quadric surface for point "
            << origin_ << " because the normal does not exist!"
            << "\nThis indicates that the input"
            << " surface mesh is of poor quality" << exit(FatalError);

        normal_ = vector::zero;
    }
    else
    {
        normal_ /= mag(normal_);
    }

    calculateBestFit();
}

//- Construct from point and other points
template<class ListType>
inline quadricFitting::quadricFitting
(
    const point& origin,
    const ListType& otherPoints
)
:
    origin_(origin),
    normal_(),
    vecX_(),
    vecY_(),
    otherPoints_(),
    transformedPoints_(),
    coefficients_()
{
    otherPoints_.setSize(otherPoints.size());
    forAll(otherPoints, i)
        otherPoints_[i] = otherPoints[i];

    calculateNormalVector();

    calculateBestFit();
}

inline quadricFitting::~quadricFitting()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const vector& quadricFitting::normal() const
{
    return normal_;
}

//- Return min curvature
inline scalar quadricFitting::minCurvature() const
{
    return
    (
        coefficients_[0] + coefficients_[1] -
        sqrt(sqr(coefficients_[0] - coefficients_[1]) + sqr(coefficients_[2]))
    );
}

//- Return max curvature
inline scalar quadricFitting::maxCurvature() const
{
    return
    (
        coefficients_[0] + coefficients_[1] +
        sqrt(sqr(coefficients_[0] - coefficients_[1]) + sqr(coefficients_[2]))
    );
}

//- Return mean curvature
inline scalar quadricFitting::meanCurvature() const
{
    return 0.5 * (minCurvature() + maxCurvature());
}

//- Return Gaussian curvature
inline scalar quadricFitting::gaussianCurvature() const
{
    return minCurvature() * maxCurvature();
}

//- Return min curvature vector
inline vector quadricFitting::minCurvatureVector() const
{
    const scalar theta =
    0.5 * Foam::atan
    (
        coefficients_[2] /
        stabilise(coefficients_[0] - coefficients_[1], VSMALL)
    );

    return Foam::cos(theta) * vecX_ + Foam::sin(theta) * vecY_;
}

//- Return max curvature vector
inline vector quadricFitting::maxCurvatureVector() const
{
    const scalar theta =
        0.5 * Foam::atan
        (
            coefficients_[2] /
            stabilise(coefficients_[0] - coefficients_[1], VSMALL)
        );

    return Foam::sin(theta) * vecX_ + Foam::cos(theta) * vecY_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
