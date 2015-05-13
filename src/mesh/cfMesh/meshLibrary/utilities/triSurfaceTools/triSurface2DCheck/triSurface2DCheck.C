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

#include "triSurface2DCheck.H"
#include "triSurfModifier.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurface2DCheck::createCovarianceMatrix()
{
    const vectorField& fNormals = surf_.facetNormals();

    //- find the normal vector of the best-fitting plane
    covarianceMatrix_ = symmTensor::zero;

    forAll(fNormals, tI)
    {
        vector fn = fNormals[tI];
        fn /= (mag(fn) + VSMALL);

        covarianceMatrix_ += symm(fn * fn);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurface2DCheck::triSurface2DCheck(const triSurf& surface)
:
    surf_(surface),
    covarianceMatrix_(symmTensor::zero)
{
    createCovarianceMatrix();
}

triSurface2DCheck::~triSurface2DCheck()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool triSurface2DCheck::is2DSurface() const
{
    const pointField& points = surf_.points();

    const vector eigenVal = eigenValues(covarianceMatrix_);

    //- the smallest eigenvalue must be zero in case all face normals
    //- lie in a plane
    if( mag(eigenVal[0]) > SMALL )
    {
        WarningIn("bool triSurface2DCheck::is2DSurface() const")
            << "Surface mesh is in 3D space!"
            << " This may result in an invalid mesh!" << endl;

        return false;
    }

    //- calculate the plane normal as a cross prduct of the two
    //- eigenVectors spanning the plane
    const vector n
    (
        eigenVector(covarianceMatrix_, eigenVal[1]) ^
        eigenVector(covarianceMatrix_, eigenVal[2])
    );

    //- check if the plane is in the x-y plane of the coordinate system
    if( mag(n.x()) > SMALL || mag(n.y()) > SMALL )
    {
        //- this could be a 2D surface, but it is not in the x-y plane
        WarningIn("bool triSurface2DCheck::is2DSurface() const")
            << "The surface mesh IS NOT IN THE X-Y PLANE!!!!"
            << " This will result in a mesh without any cells" << endl;

        return false;
    }

    //- check if the points in the 2D surface have uniform z coordinates
    boundBox bb(points);
    forAll(points, pI)
    {
        const point& p = points[pI];

        if
        (
            mag(p.z() - bb.max().z()) > SMALL &&
            mag(p.z() - bb.min().z()) > SMALL
        )
        {
            WarningIn("bool triSurface2DCheck::is2DSurface() const")
                << "z coordinates of the 2D surface are not uniform" << endl;

            return false;
        }
    }

    Info << "Detected a 2D surface in the x-y plane" << endl;

    return true;
}

void triSurface2DCheck::createSubsets()
{
    const pointField& points = surf_.points();
    const vectorField& fNormals = surf_.facetNormals();

    //- create a subset containing faces having non-zero z coordinate
    //- of the normals
    triSurf& surf = const_cast<triSurf&>(surf_);
    const label badFacetsId = surf.addFacetSubset("badFacets");
    forAll(fNormals, triI)
    {
        vector fn = fNormals[triI];
        fn /= (mag(fn) + VSMALL);

        if( mag(fn.z()) > SMALL )
            surf.addFacetToSubset(badFacetsId, triI);
    }

    //- create a subset containing points which are not
    //- in z-min and z-max planes
    const label badPointsId = surf.addPointSubset("badPointsId");
    boundBox bb(points);
    forAll(points, pI)
    {
        const point& p = points[pI];

        if
        (
            mag(p.z() - bb.max().z()) > SMALL &&
            mag(p.z() - bb.min().z()) > SMALL
        )
        {
            surf.addPointToSubset(badPointsId, pI);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
