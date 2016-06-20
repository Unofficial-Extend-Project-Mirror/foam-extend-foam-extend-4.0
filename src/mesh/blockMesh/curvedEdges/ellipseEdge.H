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

Class
    Foam::ellipseEdge

Description
    ellipseEdge class : defines the edge of an ellipse in terms of two
    points on the ellipse, its centre and eccentricity, for example:

    @verbatim
    ellipse 8 9 (32.8182156135093 0 0) -0.341427833787956
    ellipse 9 10 (0 11.9214352419999 -32.8182156135093) 0.341427833787956
    @endverbatim

    This still defines two elipses and negative eccentricity selects the
    "other" one.

SourceFiles
    ellipseEdge.C

Authors
    Henrik Rusche, Wikki GmbH
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef ellipseEdge_H
#define ellipseEdge_H

#include "curvedEdge.H"
#include "ellipticCylindricalCS.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam::mathematicalConstant;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ellipseEdge Declaration
\*---------------------------------------------------------------------------*/

class ellipseEdge
:
    public curvedEdge
{
    // Private data

        vector p1_, p2_, centre_;
        scalar eccentricity_;

        scalar theta1_, theta2_;
        scalar mu_;
        ellipticCylindricalCS cs_;

        ellipticCylindricalCS calcCS();

public:

    //- Runtime type information
    TypeName("ellipse");


    // Constructors

        //- Construct from components
        ellipseEdge
        (
            const pointField& points,
            const label start, const label end,
            const vector&,
            const scalar
        );

        //- Construct from Istream setting pointsList
        ellipseEdge(const pointField& points, Istream&);


    // Destructor

        virtual ~ellipseEdge(){}


    // Member Functions

        //- Return the position of a point on the curve given by
        //  the parameter 0 <= lambda <= 1
        vector position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
