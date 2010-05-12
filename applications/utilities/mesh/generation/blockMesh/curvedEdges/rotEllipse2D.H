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

Class
    Foam::ellipseEdge

Description
    2D ellipse with origin at (0.0, 0.0) constructed from two points
    and eccentricity

SourceFiles
    rotEllipse2D.H

Authors
    Henrik Rusche, Wikki GmbH
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef rotEllipse1D_H
#define rotEllipse2D_H

#include "mathematicalConstants.H"
#include "BisectionRoot.H"
#include "vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam::mathematicalConstant;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ellipseEdge Declaration
\*---------------------------------------------------------------------------*/

class rotEllipse2D
{
    // Ratio of minor to major axis
    scalar aOverB_;

    // Rotation around the origin
    scalar alpha_;

    // Scaling Factor
    scalar scale_;

    class func 
    {
        scalar x1_, x2_, y1_, y2_, eFactor_;

        public:

        func(scalar x1, scalar y1, scalar x2, scalar y2, scalar ecc)
        :
            x1_(x1),
            x2_(x2),
            y1_(y1),
            y2_(y2),
            eFactor_(1.0 - sqr(ecc))
        {}

        scalar operator()(const scalar a) const
        {
            scalar s = sin(a);
            scalar c = cos(a);

            return eFactor_*sqr(x1_*c + y1_*s) + sqr(x1_*s - y1_*c)
                - eFactor_*sqr(x2_*c + y2_*s) - sqr(x2_*s - y2_*c);
        }

    };

public:

    rotEllipse2D(vector2D p1, vector2D p2, scalar ecc)
    :
        aOverB_(sqrt(1.0 - sqr(ecc))),
        scale_(1.0)
    {
        func f (p1.x(), p1.y(), p2.x(), p2.y(), ecc);

        if ( mag(f(0.0)) < 1e-5 )
        {
            alpha_ = 0.0;
        }
        else if ( mag(f(piByTwo)) < 1e-5 )
        {
            alpha_ = piByTwo;
        }
        else
        {
            BisectionRoot<func> findRoot(f, 1e-5);
            if ( ecc > 0 )
            {
                alpha_ = findRoot.root(0.0, piByTwo);
            }
            else
            {
                alpha_ = findRoot.root(piByTwo, pi);
            }
        }

        scale_ = mag(p1)/mag(point(theta(p1)));
    }

    // Return the rotation of the ellipse
    scalar alpha() const
    {
        return alpha_;
    }

    // Ratio of major to minor axis
    scalar aOverB() const
    {
        return aOverB_;
    }

    // Return the scaling factor
    scalar scale() const
    {
        return scale_;
    }

    // Return the point for a given parameter
    vector2D point(const scalar theta) const
    {
        scalar xl = scale_*cos(theta);
        scalar yl = scale_*aOverB_*sin(theta);

        scalar s = sin(alpha_);
        scalar c = cos(alpha_);

        return vector2D(xl*c - yl*s, yl*c + xl*s);
    }

    // Return the parameter of the ellipse for a given point
    scalar theta(const vector2D& p) const
    {
        scalar s = sin(alpha_);
        scalar c = cos(alpha_);

        return atan2(p.y()*c - p.x()*s, aOverB_*(p.x()*c + p.y()*s));
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
