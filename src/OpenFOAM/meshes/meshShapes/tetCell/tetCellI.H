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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
inline tetCell::tetCell()
{}


//- Construct from components
inline tetCell::tetCell
(
    const label a,
    const label b,
    const label c,
    const label d
)
{
    operator[](0) = a;
    operator[](1) = b;
    operator[](2) = c;
    operator[](3) = d;
}


inline tetCell::tetCell(Istream& is)
:
    FixedList<label, 4>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline triFace tetCell::face(const label facei) const
{
    // Warning. Ordering of faces needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell
    static const label a[] = {1, 0, 0, 0};
    static const label b[] = {2, 3, 1, 2};
    static const label c[] = {3, 2, 3, 1};

#   ifdef FULLDEBUG
    if (facei >= 4)
    {
        FatalErrorIn("tetCell::tetEdge(const label facei) const")
            << "index out of range 0 -> 3. facei = " << facei
            << abort(FatalError);
    }
#   endif

    return triFace
    (
        operator[](a[facei]),
        operator[](b[facei]),
        operator[](c[facei])
    );
}


inline label tetCell::edgeFace(const label edgei) const
{
    // Warning. Ordering of faces needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell
    //static const label edgeFaces[6] = {2, 1, 1, 0, 0, 0};
    static const label edgeFaces[6] = {2, 3, 1, 0, 0, 1};

#   ifdef FULLDEBUG
    if (edgei >= 6)
    {
        FatalErrorIn
        (
            "tetCell::edgeFace(const label edgei, const label facei)"
            "const"
        )   << "edge index out of range 0 -> 5. edgei = " << edgei
            << abort(FatalError);
    }
#   endif

    return edgeFaces[edgei];
}


inline label tetCell::edgeAdjacentFace
(
    const label edgei,
    const label facei
) const
{
    // Warning. Ordering of faces needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell
    static const label adjacentFace[6][4] =
    {
        {-1, -1, 3, 2},
        {-1, 3, -1, 1},
        {-1, 2, 1, -1},
        {2, -1, 0, -1},
        {3, -1, -1, 0},
        {1, 0, -1, -1}
    };
    
#   ifdef FULLDEBUG
    if (facei >= 4)
    {
        FatalErrorIn
        (
            "tetCell::edgeAdjacentFace(const label edgei, const label facei)"
            "const"
        )   << "face index out of range 0 -> 3. facei = " << facei
            << abort(FatalError);
    }

    if (edgei >= 6)
    {
        FatalErrorIn
        (
            "tetCell::edgeAdjacentFace(const label edgei, const label facei)"
            "const"
        )   << "edge index out of range 0 -> 5. edgei = " << edgei
            << abort(FatalError);
    }
#   endif

    return adjacentFace[edgei][facei];
}


inline edge tetCell::tetEdge(const label edgei) const
{
    // Warning. Ordering of edges needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell
    // 
    static const label start[] = {0, 0, 0, 3, 1, 3};
    static const label end[] = {1, 2, 3, 1, 2, 2};

#   ifdef FULLDEBUG
    if (edgei >= 6)
    {
        FatalErrorIn("tetCell::tetEdge(const label edgei) const")
            << "index out of range 0 -> 5. edgei = " << edgei
            << abort(FatalError);
    }
#   endif

    return edge(operator[](start[edgei]), operator[](end[edgei]));
}


inline tetPointRef tetCell::tet(const pointField& points) const
{
    return tetPointRef
    (
        points[operator[](0)],
        points[operator[](1)],
        points[operator[](2)],
        points[operator[](3)]
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
