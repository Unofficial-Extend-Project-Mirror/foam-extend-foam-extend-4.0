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

\*---------------------------------------------------------------------------*/

#include "boundBox.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::boundBox::great(VGREAT);

const Foam::boundBox Foam::boundBox::greatBox
(
    point(-VGREAT, -VGREAT, -VGREAT),
    point(VGREAT, VGREAT, VGREAT)
);


const Foam::boundBox Foam::boundBox::invertedBox
(
    point(VGREAT, VGREAT, VGREAT),
    point(-VGREAT, -VGREAT, -VGREAT)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundBox::boundBox(const pointField& points, const bool doReduce)
:
    min_(vector::zero),
    max_(vector::zero)
{
    if (points.size() == 0)
    {
        if (Pstream::parRun() && doReduce)
        {
            // Use values which get overwritten by reduce minOp,maxOp below
            min_ = point(VGREAT, VGREAT, VGREAT);
            max_ = point(-VGREAT, -VGREAT, -VGREAT);
        }
        else
        {
            WarningIn("boundBox::boundBox(const pointField& points)")
                << "Cannot find bounding box for zero sized pointField, "
                   "returning zero"
                << endl;

            return;
        }
    }
    else
    {
        min_ = points[0];
        max_ = points[0];

        forAll(points, i)
        {
            min_ = ::Foam::min(min_, points[i]);
            max_ = ::Foam::max(max_, points[i]);
        }
    }

    if (doReduce)
    {
        // Reduce parallel information
        reduce(min_, minOp<point>());
        reduce(max_, maxOp<point>());
    }
}


boundBox::boundBox(Istream& is)
:
    min_(is),
    max_(is)
{}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const boundBox& b)
{
    return os << b.min() << token::SPACE << b.max();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
