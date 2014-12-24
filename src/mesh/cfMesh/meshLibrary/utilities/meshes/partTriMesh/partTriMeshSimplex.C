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

#include "Map.H"
#include "partTriMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTriMeshSimplex::partTriMeshSimplex
(
    const partTriMesh& tm,
    const label pI
)
:
    pts_(),
    trias_()
{
    const pointField& points = tm.points();
    const LongList<labelledTri>& trias = tm.triangles();
    const VRWGraph& pt = tm.pointTriangles();
    
    trias_.setSize(pt.sizeOfRow(pI));
    label counter(0);

    Map<label> addr(2*pt.sizeOfRow(pI));
    forAllRow(pt, pI, tI)
    {
        const labelledTri& tri = trias[pt(pI, tI)];
        for(label i=0;i<3;++i)
        {
            const label tpI = tri[i];

            if( !addr.found(tpI) )
            {
                addr.insert(tpI, counter);
                pts_.append(points[tpI]);
                ++counter;
            }
        }
        
        # ifdef DEBUGSmooth
        Info << "Tet " << tetI << " is " << tet << endl;
        # endif
        
        label pos(-1);
        for(label i=0;i<3;++i)
            if( tri[i] == pI )
            {
                pos = i;
                break;
            }

        switch( pos )
        {
            case 0:
            {
                trias_[tI] = triFace(addr[tri[0]], addr[tri[1]], addr[tri[2]]);
            } break;
            case 1:
            {
                trias_[tI] = triFace(addr[tri[1]], addr[tri[2]], addr[tri[0]]);
            } break;
            case 2:
            {
                trias_[tI] = triFace(addr[tri[2]], addr[tri[0]], addr[tri[1]]);
            } break;
            default:
            {
                FatalErrorIn
                (
                    "partTriMeshSimplex::partTriMeshSimplex("
                    "(const partTriMesh& tm, const label pI)"
                ) << "Point " << pI << " is not present in triangle" << tri
                    << abort(FatalError);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
partTriMeshSimplex::~partTriMeshSimplex()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vector partTriMeshSimplex::normal() const
{
    vector normal(vector::zero);
    scalar magN(0.0);

    forAll(trias_, tI)
    {
        const triFace& t = trias_[tI];

        vector n
        (
            0.5 * ((pts_[t[1]] - pts_[t[0]]) ^ (pts_[t[2]] - pts_[t[0]]))
        );
        const scalar magn = mag(n);

        normal += n;
        magN += magn;
    }

    return (normal / (magN + VSMALL));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
