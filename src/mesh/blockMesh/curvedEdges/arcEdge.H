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
    Foam::arcEdge

Description
    Defines the arcEdge of a circle in terms of 3 points on its circumference

SourceFiles
    arcEdge.C

\*---------------------------------------------------------------------------*/

#ifndef arcEdge_H
#define arcEdge_H

#include "curvedEdge.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class arcEdge Declaration
\*---------------------------------------------------------------------------*/

class arcEdge
:
    public curvedEdge
{
    // Private data

        point p1_, p2_, p3_;
        cylindricalCS cs_;
        scalar angle_;
        scalar radius_;

    // Private Member Functions

        //- Calculate the coordinate system, angle and radius
        cylindricalCS calcAngle();

        //- Disallow default bitwise copy construct
        arcEdge(const arcEdge&);

        //- Disallow default bitwise assignment
        void operator=(const arcEdge&);

public:

    //- Runtime type information
    TypeName("arc");


    // Constructors

        //- Construct from components
        arcEdge
        (
            const pointField& points,
            const label start, const label end,
            const point& pMid
        );

        //- Construct from Istream setting pointsList
        arcEdge(const pointField& points, Istream&);


    //- Destructor
    virtual ~arcEdge(){}


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
