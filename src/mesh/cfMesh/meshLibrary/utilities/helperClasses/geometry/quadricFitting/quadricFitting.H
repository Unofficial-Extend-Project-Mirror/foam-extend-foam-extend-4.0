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

Class
    quadricFitting

Description
    A class for fitting a quadric surface through the set of points. It is
    useful for curvature estimation and for calculating normals.

SourceFiles
    quadricFittingI.H

\*---------------------------------------------------------------------------*/

#ifndef quadricFitting_H
#define quadricFitting_H

#include "DynList.H"
#include "plane.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class quadricFitting Declaration
\*---------------------------------------------------------------------------*/

class quadricFitting
{
    // Private data
        //- seed point for fitting the quadric surface
        point origin_;

        //- normal vector at the origin
        vector normal_;

        //- x-coordinate vector for the transformed coordinate system
        vector vecX_;

        //- y-coordinate vector for the transformed coordinate system
        vector vecY_;

        //- other points
        DynList<point> otherPoints_;

        //- transformed coordinates of other points
        DynList<point> transformedPoints_;

        //- coefficients of the quadric surface
        scalarField coefficients_;

    // Private member functions
        //- initial estimate of the normal vector
        void calculateNormalVector();

        //- calculate transformed coordinate system
        void calculateCoordinateSystem();

        //- calculate coefficients of the quadric surface fit
        void calculateQuadricCoeffs();

        //- iterative procedure for updating normals
        //- until the best fit is achieved
        void calculateBestFit();

public:

    // Constructors

        //- Construct from point, normal, and neighbouring points
        template<class ListType>
        inline quadricFitting(const point&, const vector, const ListType&);

        //- Construct from point and other points
        template<class ListType>
        inline quadricFitting(const point&, const ListType&);

    // Destructor
        ~quadricFitting();


    // Member Functions

        // Access

            //- Return surface normal
            inline const vector& normal() const;

            //- Return min curvature
            inline scalar minCurvature() const;

            //- Return max curvature
            inline scalar maxCurvature() const;

            //- Return mean curvature
            inline scalar meanCurvature() const;

            //- Return Gaussian curvature
            inline scalar gaussianCurvature() const;

            //- Return min curvature vector
            inline vector minCurvatureVector() const;

            //- Return max curvature vector
            inline vector maxCurvatureVector() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "quadricFittingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
