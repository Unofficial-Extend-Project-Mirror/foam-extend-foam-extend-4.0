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
#include "partTetMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTetMeshSimplex::partTetMeshSimplex
(
    const partTetMesh& tm,
    const label pI
)
:
    pts_(),
    tets_()
{
    const LongList<point>& points = tm.points();
    const LongList<partTet>& tets = tm.tets();
    const VRWGraph& pt = tm.pointTets();
    
    tets_.setSize(pt.sizeOfRow(pI));
    label counter(0);

    Map<label> addr(2*pt.sizeOfRow(pI));
    forAllRow(pt, pI, tetI)
    {
        const partTet& tet = tets[pt(pI, tetI)];
        for(label i=0;i<4;++i)
        {
            const label tpI = tet[i];
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
        
        const label pos = tet.whichPosition(pI);
        switch( pos )
        {
            case 0:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet.b()],
                        addr[tet.d()],
                        addr[tet.c()],
                        addr[tet.a()]
                    );
            } break;
            case 1:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet.a()],
                        addr[tet.c()],
                        addr[tet.d()],
                        addr[tet.b()]
                    );
            } break;
            case 2:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet.a()],
                        addr[tet.d()],
                        addr[tet.b()],
                        addr[tet.c()]
                    );
            } break;
            case 3:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet.a()],
                        addr[tet.b()],
                        addr[tet.c()],
                        addr[tet.d()]
                    );
            } break;
            default:
            {
                FatalErrorIn
                (
                    "partTetMeshSimplex::partTetMeshSimplex("
                    "(const partTetMesh& tm, const label pI)"
                ) << "Point " << pI << " is not present in tet" << tet
                    << abort(FatalError);
            }
        }
    }
}

partTetMeshSimplex::partTetMeshSimplex
(
    const DynList<parPartTet>& pt,
    const label gpI
)
:
    pts_(),
    tets_()
{
    tets_.setSize(pt.size());
    label pI(0);
    
    Map<label> addr;
    forAll(pt, tetI)
    {
        const parPartTet& tet = pt[tetI];
        
        label pos(-1);
        for(label i=0;i<4;++i)
        {
            if( !addr.found(tet[i].pointLabel()) )
            {
                addr.insert(tet[i].pointLabel(), pI);
                pts_.append(tet[i].coordinates());
                ++pI;
            }
            
            if( tet[i].pointLabel() == gpI )
                pos = i;
        }
        
        switch( pos )
        {
            case 0:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet[1].pointLabel()],
                        addr[tet[3].pointLabel()],
                        addr[tet[2].pointLabel()],
                        addr[tet[0].pointLabel()]
                    );
            } break;
            case 1:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet[0].pointLabel()],
                        addr[tet[2].pointLabel()],
                        addr[tet[3].pointLabel()],
                        addr[tet[1].pointLabel()]
                    );
            } break;
            case 2:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet[0].pointLabel()],
                        addr[tet[3].pointLabel()],
                        addr[tet[1].pointLabel()],
                        addr[tet[2].pointLabel()]
                    );
            } break;
            case 3:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet[0].pointLabel()],
                        addr[tet[1].pointLabel()],
                        addr[tet[2].pointLabel()],
                        addr[tet[3].pointLabel()]
                    );
            } break;
            default:
            {
                FatalErrorIn
                (
                    "partTetMeshSimplex::partTetMeshSimplex("
                    "(const partTetMesh& tm, const label pI)"
                ) << "Point " << gpI << " is not present in tet" << tet
                    << abort(FatalError);
            }
        }
    }
}

partTetMeshSimplex::partTetMeshSimplex
(
    const DynList<point, 128>& pts,
    const DynList<partTet, 128>& tets,
    const label pointI
)
:
    pts_(pts),
    tets_(tets.size())
{
    forAll(tets, tetI)
    {
        const partTet& tet = tets[tetI];

        const label pos = tet.whichPosition(pointI);

        switch( pos )
        {
            case 0:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.b(),
                        tet.d(),
                        tet.c(),
                        tet.a()
                    );
            } break;
            case 1:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.c(),
                        tet.d(),
                        tet.b()
                    );
            } break;
            case 2:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.d(),
                        tet.b(),
                        tet.c()
                    );
            } break;
            case 3:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.b(),
                        tet.c(),
                        tet.d()
                    );
            } break;
            default:
            {
                FatalErrorIn
                (
                    "partTetMeshSimplex::partTetMeshSimplex"
                    "(const DynList<point, 128>& pts,"
                    "const DynList<partTet, 128>& tets, const label pointI)"
                ) << "Point " << pointI << " is not present in tet" << tet
                    << abort(FatalError);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
partTetMeshSimplex::~partTetMeshSimplex()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
