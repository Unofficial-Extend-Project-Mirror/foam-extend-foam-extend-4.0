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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label boundaryLayers::findNewNodeLabel
(
    const label pointI,
    const label pKey
) const
{
    if( otherVrts_.find(pointI) != otherVrts_.end() )
    {
        const std::map
        <
            label, std::map<std::pair<label, label>, label>
        >::const_iterator it = otherVrts_.find(pointI);
        const std::map<std::pair<label, label>, label>& m = it->second;

        std::map<std::pair<label, label>, label>::const_iterator mit;
        if( m.size() == 2 )
        {
            for(mit=m.begin();mit!=m.end();++mit)
            {
                if( mit->first.first != pKey )
                    return mit->second;
            }
        }
        else
        {
            for(mit=m.begin();mit!=m.end();++mit)
            {
                if( mit->first.first == pKey )
                    continue;
                if( mit->first.second == pKey )
                    continue;
                if( mit->first.first == mit->first.second )
                    continue;
                
                return mit->second;
            }
        }
    }

    return newLabelForVertex_[pointI];
}

inline void boundaryLayers::createNewCellFromEdge
(
    const edge& e,
    const label pKeyI,
    const label pKeyJ,
    FixedList<FixedList<label, 4>, 6>& cellFaces
) const
{
    const std::map<std::pair<label, label>, label>& ms =
        otherVrts_.find(e.start())->second;
    const std::map<std::pair<label, label>, label>& me =
        otherVrts_.find(e.end())->second;
    
    # ifdef DEBUGLayer
    Info << "Creating cell for edge " << edgeI << " with nodes " << e << endl;
    Info << "pKeyI " << pKeyI << endl;
    Info << "pKeyJ " << pKeyJ << endl;
    std::map<std::pair<label, label>, label>::const_iterator iter;
    for(iter=ms.begin();iter!=ms.end();++iter)
        Info << "1. Pair (" << iter->first.first << ", " 
            << iter->first.second << ") has value " << iter->second << endl;
    for(iter=me.begin();iter!=me.end();++iter)
        Info << "2. Pair (" << iter->first.first << ", " 
            << iter->first.second << ") has value " << iter->second << endl;
    # endif
    
    label p0s(-1), p1s(-1), ns(-1), p0e(-1), p1e(-1), ne(-1);
    
    if( ms.size() == 2 )
    {
        p0s = ms.find(std::pair<label, label>(pKeyI, pKeyI))->second;
        p1s = ms.find(std::pair<label, label>(pKeyJ, pKeyJ))->second;
        ns = newLabelForVertex_[e.start()];
    }
    else
    {
        std::map<std::pair<label, label>, label>::const_iterator it;
        for(it=ms.begin();it!=ms.end();++it)
        {
            if(
                (it->first.first != pKeyJ) && (it->first.second != pKeyJ)
                && (it->first.first != it->first.second)
            )
            {
                p0s = it->second;
            }
            else if(
                (it->first.first != pKeyI) && (it->first.second != pKeyI)
                && (it->first.first != it->first.second)
            )
            {
                p1s = it->second;
            }
            else if(
                (it->first.first == it->first.second) &&
                (it->first.first != pKeyI) && (it->first.first != pKeyJ)
            )
            {
                ns = it->second;
            }
        }
    }
    if( me.size() == 2 )
    {
        p0e = me.find(std::pair<label, label>(pKeyI, pKeyI))->second;
        p1e = me.find(std::pair<label, label>(pKeyJ, pKeyJ))->second;
        ne = newLabelForVertex_[e.end()];
    }
    else
    {
        std::map<std::pair<label, label>, label>::const_iterator it;
        for(it=me.begin();it!=me.end();++it)
        {
            if(
                (it->first.first != pKeyJ) && (it->first.second != pKeyJ)
                && (it->first.first != it->first.second)
            )
            {
                p0e = it->second;
            }
            else if(
                (it->first.first != pKeyI) && (it->first.second != pKeyI)
                && (it->first.first != it->first.second)
            )
            {
                p1e = it->second;
            }
            else if(
                (it->first.first == it->first.second) &&
                (it->first.first != pKeyI) && (it->first.first != pKeyJ)
            )
            {
                ne = it->second;
            }
        }
    }
    
    //- F0
    cellFaces[0][0] = e.end();
    cellFaces[0][1] = e.start();
    cellFaces[0][2] = p1s;
    cellFaces[0][3] = p1e;

    //- F1
    cellFaces[1][0] = p0e;
    cellFaces[1][1] = ne;
    cellFaces[1][2] = ns;
    cellFaces[1][3] = p0s;

    //- F2
    cellFaces[2][0] = e.start();
    cellFaces[2][1] = e.end();
    cellFaces[2][2] = p0e;
    cellFaces[2][3] = p0s;

    //- F3
    cellFaces[3][0] = p1s;
    cellFaces[3][1] = ns;
    cellFaces[3][2] = ne;
    cellFaces[3][3] = p1e;

    //- F4
    cellFaces[4][0] = e.end();
    cellFaces[4][1] = p1e;
    cellFaces[4][2] = ne;
    cellFaces[4][3] = p0e;

    //- F5
    cellFaces[5][0] = e.start();
    cellFaces[5][1] = p0s;
    cellFaces[5][2] = ns;
    cellFaces[5][3] = p1s;
}

inline void boundaryLayers::createNewCellFromNode
(
    const label pointI,
    const DynList<label, 3>& pKeys,
    FixedList<FixedList<label, 4>, 6>& cellFaces
) const
{
    const std::map<std::pair<label, label>, label>& m =
        otherVrts_.find(pointI)->second;
    
    //- create labels before creating cells
    const label n = newLabelForVertex_[pointI];
    const label p00 =
        m.find(std::pair<label, label>(pKeys[0], pKeys[0]))->second;
    const label p11 =
        m.find(std::pair<label, label>(pKeys[1], pKeys[1]))->second;
    const label p22 =
        m.find(std::pair<label, label>(pKeys[2], pKeys[2]))->second;

    std::pair<label, label> pr;
    pr.first = pKeys[0];
    pr.second = pKeys[1];
    if( m.find(pr) == m.end() )
    {
        pr.first = pKeys[1];
        pr.second = pKeys[0];
    }
    const label p01 = m.find(pr)->second;
    
    pr.first = pKeys[0];
    pr.second = pKeys[2];
    if( m.find(pr) == m.end() )
    {
        pr.first = pKeys[2];
        pr.second = pKeys[0];
    }
    const label p02 = m.find(pr)->second;
    
    pr.first = pKeys[1];
    pr.second = pKeys[2];
    if( m.find(pr) == m.end() )
    {
        pr.first = pKeys[2];
        pr.second = pKeys[1];
    }
    const label p12 = m.find(pr)->second;
    
    //- create the cell and append it
    
    //- F0
    cellFaces[0][0] = pointI;
    cellFaces[0][1] = p02;
    cellFaces[0][2] = p00;
    cellFaces[0][3] = p01;
    //- F1
    cellFaces[1][0] = p12;
    cellFaces[1][1] = p11;
    cellFaces[1][2] = n;
    cellFaces[1][3] = p22;

    //- F2
    cellFaces[2][0] = pointI;
    cellFaces[2][1] = p01;
    cellFaces[2][2] = p11;
    cellFaces[2][3] = p12;
    //- F3
    cellFaces[3][0] = p02;
    cellFaces[3][1] = p22;
    cellFaces[3][2] = n;
    cellFaces[3][3] = p00;
    //- F4
    cellFaces[4][0] = pointI;
    cellFaces[4][1] = p12;
    cellFaces[4][2] = p22;
    cellFaces[4][3] = p02;
    //- F5
    cellFaces[5][0] = p01;
    cellFaces[5][1] = p00;
    cellFaces[5][2] = n;
    cellFaces[5][3] = p11;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
