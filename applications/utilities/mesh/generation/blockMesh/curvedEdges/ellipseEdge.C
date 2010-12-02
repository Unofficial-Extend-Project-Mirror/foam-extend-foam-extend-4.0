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

Description
    ellipseEdge class : defines the edge of an ellipse in terms of its
    centre, 2 point and eccentricity

\*---------------------------------------------------------------------------*/

#include "ellipseEdge.H"
#include "mathematicalConstants.H"
#include "EulerCoordinateRotation.H"
#include "coordinateSystem.H"
#include "curvedEdges/rotEllipse2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ellipseEdge, 0);

    // Add the curvedEdge constructor functions to the hash tables
    curvedEdge::addIstreamConstructorToTable<ellipseEdge>
        addEllipseEdgeIstreamConstructorToTable_;
}

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::ellipticCylindricalCS Foam::ellipseEdge::calcCS()
{
    if ( debug )
    {
        Info << "p1 = " << p1_ << " p2 = " << p2_ << endl;
    }

    vector p = p1_ - centre_;
    vector q = p2_ - centre_;

    if ( debug )
    {
        Info << "mag(p) = " << mag(p) << " mag(q) = " << mag(q) << endl;
    }

    vector normal = (p ^ q);
    normal /= mag(normal);

    scalar scale = mag(p);
    p /= scale;
    q /= scale;

    if ( debug )
    {
        Info << "p = " << p << " q = " << q << " n = " << normal << endl;
    }

    tensor A(p, normal ^ p, normal);

    vector2D p2D(1.0, 0.0);
    vector qTmp = (A & q);
    vector2D q2D(qTmp.x(), qTmp.y());

    if ( debug )
    {
        Info << "p2D = " << p2D << endl;
        Info << "q2D = " << q2D << endl;
    }

    rotEllipse2D ellipse(p2D, q2D, eccentricity_);

    scalar alpha = ellipse.alpha();
    theta1_ = ellipse.theta(p2D);
    theta2_ = ellipse.theta(q2D);
    mu_ = atanh(ellipse.aOverB());

    scalar a = scale*ellipse.scale();

    if ( debug )
    {
        Info << "alpha    = " << alpha << endl;
        Info << "scale    = " << ellipse.scale() << endl;
        Info << "a        = " << a << endl;
        Info << "theta1   = " << theta1_ << " theta2 = " << theta2_ << endl;

        vector2D pn = ellipse.point(theta1_);
        vector2D qn = ellipse.point(theta2_);

        Info << "pn = " << pn << tab << pn/mag(pn) << " mag = " << mag(pn) << endl;
        Info << "qn = " << qn << tab << qn/mag(pn) << " mag = " << mag(qn) << endl;

        Info << (pn & qn)/mag(pn)/mag(qn) << tab << atan2(qn.y(), qn.x()) << endl;
    }

    // set up and return the local coordinate system
    return ellipticCylindricalCS
    (
        "tmpCS",
        centre_,
        (vector(0,0,1) & A),
        (vector(cos(alpha), sin(alpha), 0) & A),
        a/cosh(mu_)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::ellipseEdge::ellipseEdge
(
    const pointField& points,
    const label start,
    const label end,
    const vector& centre,
    const scalar eccentricity
)
:
    curvedEdge(points, start, end),
    p1_(points_[start_]),
    p2_(points_[end_]),
    centre_(centre),
    eccentricity_(eccentricity),
    cs_(calcCS())
{}


// Construct from Istream
Foam::ellipseEdge::ellipseEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    p1_(points_[start_]),
    p2_(points_[end_]),
    centre_(is),
    eccentricity_(readScalar(is)),
    cs_(calcCS())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::ellipseEdge::position(const scalar lambda) const
{
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("ellipseEdge::position(const scalar lambda) const")
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    if (lambda < SMALL)
    {
        return p1_;
    }
    else if (lambda > 1-SMALL)
    {
        return p2_;
    }
    else
    {
        scalar theta = theta1_ + lambda*(theta2_-theta1_);
        return cs_.globalPosition(vector(mu_, theta/pi*180, 0.0));
    }
}


//- Return the length of the curve
Foam::scalar Foam::ellipseEdge::length() const
{
    notImplemented("ellipseEdge::length() const");
    return 1.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
