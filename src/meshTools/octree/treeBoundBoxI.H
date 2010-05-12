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

#include "treeBoundBox.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null setting points to zero
inline treeBoundBox::treeBoundBox()
:
    boundBox()
{}


// Construct from components
inline treeBoundBox::treeBoundBox(const point& min, const point& max)
:
    boundBox(min, max)
{}


// Construct from components
inline treeBoundBox::treeBoundBox(const boundBox& bb)
:
    boundBox(bb)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar treeBoundBox::minDim() const
{
    return ::Foam::min
    (
        max().x() - min().x(),
        ::Foam::min
        (
            max().y() - min().y(),
            max().z() - min().z()
        )
    );
}


inline scalar treeBoundBox::maxDim() const
{
    return ::Foam::max
    (
        max().x() - min().x(),
        ::Foam::max
        (
            max().y() - min().y(),
            max().z() - min().z()
        )
    );
}


inline scalar treeBoundBox::avgDim() const
{
    return
    (
        (max().x() - min().x()) +
        (max().y() - min().y()) +
        (max().z() - min().z())
    )/3.0;
}


inline scalar treeBoundBox::typDim() const
{
    return avgDim();
}


inline point treeBoundBox::mid() const
{
    return 0.5*(min() + max());
}


inline point treeBoundBox::corner(const direction octant) const
{
    return point
    (
        (octant&RIGHTHALF) ? max().x() : min().x(),
        (octant&TOPHALF)   ? max().y() : min().y(),
        (octant&FRONTHALF) ? max().z() : min().z()
    );
}

// Returns octant in which sample resides. Reverse of subBbox.
inline direction treeBoundBox::subOctant(const point& sample) const
{
    point mid = 0.5*(max() + min());

    direction octant = 0;

    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }

    return octant;
}


// Returns octant in which sample resides. Reverse of subBbox. Precalculated
// midpoint
inline direction treeBoundBox::subOctant
(
    const point& mid,
    const point& sample
)
{
    direction octant = 0;

    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }

    return octant;
}


// Returns octant in which sample resides. Reverse of subBbox. Flags sample
// exactly on edge.
inline direction treeBoundBox::subOctant(const point& sample, bool& onEdge)
 const
{
    point mid = 0.5*(max() + min());

    direction octant = 0;

    onEdge = false;
    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }
    else if (sample.x() == mid.x())
    {
        onEdge = true;
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }
    else if (sample.y() == mid.y())
    {
        onEdge = true;
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }
    else if (sample.z() == mid.z())
    {
        onEdge = true;
    }

    return octant;
}


// Returns octant in which sample resides. Reverse of subBbox. Precalculated
// midpoint
inline direction treeBoundBox::subOctant
(
    const point& mid,
    const point& sample,
    bool& onEdge
)
{
    direction octant = 0;

    onEdge = false;
    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }
    else if (sample.x() == mid.x())
    {
        onEdge = true;
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }
    else if (sample.y() == mid.y())
    {
        onEdge = true;
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }
    else if (sample.z() == mid.z())
    {
        onEdge = true;
    }

    return octant;
}


// Returns octant in which intersection resides.
// Precalculated midpoint. If the sample is on the dividing line between
// the octants the direction vector determines which octant to use
// (i.e. in which octant the sample would be if it were moved along dir)
inline direction treeBoundBox::subOctant
(
    const point& mid,
    const vector& dir,
    const point& sample,
    bool& onEdge
)
{
    direction octant = 0;

    onEdge = false;

    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }
    else if (sample.x() == mid.x())
    {
        onEdge = true;

        if (dir.x() > 0)
        {
            octant |= treeBoundBox::RIGHTHALF;
        }
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }
    else if (sample.y() == mid.y())
    {
        onEdge = true;

        if (dir.y() > 0)
        {
            octant |= treeBoundBox::TOPHALF;
        }
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }
    else if (sample.z() == mid.z())
    {
        onEdge = true;

        if (dir.z() > 0)
        {
            octant |= treeBoundBox::FRONTHALF;
        }
    }

    return octant;
}


// Returns reference to octantOrder which defines the 
// order to do the search.
inline void treeBoundBox::searchOrder
(
    const point& sample,
    FixedList<direction,8>& octantOrder
) const
{
    vector dist = mid() - sample;

    direction octant = 0;

    if (dist.x() < 0)
    {
        octant |= treeBoundBox::RIGHTHALF;
        dist.x() *= -1;
    }

    if (dist.y() < 0)
    {
        octant |= treeBoundBox::TOPHALF;
        dist.y() *= -1;
    }

    if (dist.z() < 0)
    {
        octant |= treeBoundBox::FRONTHALF;
	dist.z() *= -1;
    }

    direction min = 0;
    direction mid = 0;
    direction max = 0;

    if( dist.x() < dist.y())
    {
        if( dist.y() < dist.z())
	{
	    min = treeBoundBox::RIGHTHALF;
	    mid = treeBoundBox::TOPHALF;
	    max = treeBoundBox::FRONTHALF;
	}
	else if( dist.z() < dist.x())
	{
	    min = treeBoundBox::FRONTHALF;
            mid = treeBoundBox::RIGHTHALF;
	    max = treeBoundBox::TOPHALF;
	}
	else
	{
	    min = treeBoundBox::RIGHTHALF;
	    mid = treeBoundBox::FRONTHALF;
	    max = treeBoundBox::TOPHALF;
	}
    }
    else
    {
        if( dist.z() < dist.y())
	{
	    min = treeBoundBox::FRONTHALF;
	    mid = treeBoundBox::TOPHALF;
	    max = treeBoundBox::RIGHTHALF;	    
	}
	else if( dist.x() < dist.z())
	{
	    min = treeBoundBox::TOPHALF;
	    mid = treeBoundBox::RIGHTHALF;
	    max = treeBoundBox::FRONTHALF;
	}
	else
	{
	    min = treeBoundBox::TOPHALF;
	    mid = treeBoundBox::FRONTHALF;
	    max = treeBoundBox::RIGHTHALF;
	}
    }
    // Primary subOctant
    octantOrder[0] = octant;
    // subOctants joined to the primary by faces.
    octantOrder[1] = octant ^ min;
    octantOrder[2] = octant ^ mid;
    octantOrder[3] = octant ^ max;
    // subOctants joined to the primary by edges.
    octantOrder[4] = octantOrder[1] ^ mid;
    octantOrder[5] = octantOrder[1] ^ max;
    octantOrder[6] = octantOrder[2] ^ max;
    // subOctants joined to the primary by corner.
    octantOrder[7] = octantOrder[4] ^ max;
}


inline bool treeBoundBox::intersects(const treeBoundBox& bb) const
{
    return boundBox::intersects(bb);
}

// true if bb's intersect or overlap.
// Note: <= to make sure we catch all.
inline bool treeBoundBox::overlaps(const treeBoundBox& bb) const
{
    return boundBox::overlaps(bb);
}

inline bool treeBoundBox::contains(const point& sample) const
{
    return
    (
        (sample.x() >= min().x()) &&
        (sample.y() >= min().y()) &&
        (sample.z() >= min().z()) &&
        (sample.x() <= max().x()) &&
        (sample.y() <= max().y()) &&
        (sample.z() <= max().z())
    );
}


//- Return slightly wider bounding box
inline treeBoundBox treeBoundBox::extend(Random& rndGen, const scalar s) const
{
    treeBoundBox bb(*this);

    vector span(bb.max() - bb.min());

    // Make 3D
    scalar magSpan = Foam::mag(span);

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        span[dir] = Foam::max(s*magSpan, span[dir]);
    }

    bb.min() -= cmptMultiply(s*rndGen.vector01(), span);
    bb.max() += cmptMultiply(s*rndGen.vector01(), span);

    return bb;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
