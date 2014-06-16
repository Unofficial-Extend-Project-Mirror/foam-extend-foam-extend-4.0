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

#include "triangle.H"
#include "tetrahedron.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline partTet::partTet()
{
}
    
inline partTet::partTet
(
    const label a,
    const label b,
    const label c,
    const label d
)
{
    data_[0] = a;
    data_[1] = b;
    data_[2] = c;
    data_[3] = d;
}

inline partTet::~partTet()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline label partTet::a() const
{
    return data_[0];
}

inline label partTet::b() const
{
    return data_[1];
}

inline label partTet::c() const
{
    return data_[2];
}

inline label partTet::d() const
{
    return data_[3];
}

inline label partTet::size() const
{
    return 4;
}

inline label partTet::whichPosition(const label pointI) const
{
    for(label i=0;i<4;++i)
        if( data_[i] == pointI )
            return i;
        
    return -1;
}

template<class PointField>
inline vector partTet::Sa(const PointField& points) const
{
    triangle<point, point> tria
    (
        points[data_[1]],
        points[data_[2]],
        points[data_[3]]
    );
    
    return tria.normal();
    //return triangle<point, point>(b_, c_, d_).normal();
}

template<class PointField>
inline vector partTet::Sb(const PointField& points) const
{
    triangle<point, point> tria
    (
        points[data_[0]],
        points[data_[3]],
        points[data_[2]]
    );
    
    return tria.normal();
    //return triangle<point, point>(a_, d_, c_).normal();
}

template<class PointField>
inline vector partTet::Sc(const PointField& points) const
{
    triangle<point, point> tria
    (
        points[data_[0]],
        points[data_[1]],
        points[data_[3]]
    );
    
    return tria.normal();
    //return triangle<point, point>(a_, b_, d_).normal();
}

template<class PointField>
inline vector partTet::Sd(const PointField& points) const
{
    triangle<point, point> tria
    (
        points[data_[0]],
        points[data_[2]],
        points[data_[1]]
    );
    
    return tria.normal();
    //return triangle<point, point>(a_, c_, b_).normal();
}

template<class PointField>
inline scalar partTet::mag(const PointField& points) const
{
    tetrahedron<point, point> tet
    (
        points[data_[0]],
        points[data_[1]],
        points[data_[2]],
        points[data_[3]]
    );
    
    return tet.mag();
    //return (1.0/6.0)*(((b_ - a_) ^ (c_ - a_)) & (d_ - a_));
}

template<class PointField>
inline point partTet::crcmCentre(const PointField& points) const
{
    tetrahedron<point, point> tet
    (
        points[data_[0]],
        points[data_[1]],
        points[data_[2]],
        points[data_[3]]
    );
    
    return tet.circumCentre();
}

template<class PointField>
inline scalar partTet::crcmRadius(const PointField& points) const
{
    tetrahedron<point, point> tet
    (
        points[data_[0]],
        points[data_[1]],
        points[data_[2]],
        points[data_[3]]
    );
    
    return tet.circumRadius();
}

template<class PointField>
inline point partTet::centroid(const PointField& points) const
{
    point p = points[data_[0]];
    for(label i=1;i<4;++i)
        p += points[data_[i]];
    
    p /= 4;
    return p;
}

inline edgeList partTet::edges() const
{
    edgeList edges(6);
    edges[0] = edge(data_[0], data_[1]);
    edges[1] = edge(data_[0], data_[2]);
    edges[2] = edge(data_[0], data_[3]);
    edges[3] = edge(data_[3], data_[1]);
    edges[4] = edge(data_[1], data_[2]);
    edges[5] = edge(data_[3], data_[2]);

    return edges;
}

inline label partTet::operator[](const label i) const
{
    return data_[i];
}

inline void partTet::operator=(const partTet& tet)
{
    for(label i=0;i<4;++i)
        data_[i] = tet.data_[i];
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Ostream& operator<<(Ostream& os, const partTet& t)
{
    os  << nl
        << token::BEGIN_LIST
        << t.a() << token::SPACE << t.b()
        << token::SPACE << t.c() << token::SPACE << t.d()
        << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
