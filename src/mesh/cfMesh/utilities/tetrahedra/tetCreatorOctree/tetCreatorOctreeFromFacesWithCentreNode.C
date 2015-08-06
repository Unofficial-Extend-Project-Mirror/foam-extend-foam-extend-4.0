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

#include "tetCreatorOctree.H"
#include "demandDrivenData.H"

//#define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetCreatorOctree::createTetsFromFacesWithCentreNode()
{
    Info << "Creating tets from faces with centre node" << endl;
    
    const labelList& cubeLabel = *cubeLabelPtr_;
    const VRWGraph& nodeLabels = octreeCheck_.nodeLabels();
    const VRWGraph& subNodeLabels = *subNodeLabelsPtr_;
    const FRWGraph<label, 8>& pointLeaves = octreeCheck_.nodeLeaves();
    
    if( !faceCentreLabelPtr_ )
        faceCentreLabelPtr_ = new VRWGraph(cubeLabel.size());
    VRWGraph& faceCentreLabel = *faceCentreLabelPtr_;
    
    //- start creating tets
    forAll(pointLeaves, pointI)
    {
        label pl[8];
        bool create(true);
        
        for(label plI=0;plI<8;++plI)
        {
            pl[plI] = pointLeaves(pointI, plI);
            if( pl[plI] == -1 ) 
            {
                create = false;
                break;
            }
        }
            
        if( !create )
            continue;
        
        //- create 6 tets for each possible combination
        //- there are 12 possible combinations
        for(label fI=0;fI<6;++fI)
        {
            const label* fEdges = meshOctreeCubeCoordinates::faceEdges_[fI];
            
            for(label feI=0;feI<2;++feI)
            {
                const label feJ = feI + 2;
                
                //- the are two possible combinations of edges for each face
                const label* sEdge =
                    meshOctreeCubeCoordinates::edgeNodes_[fEdges[feI]];
                const label* eEdge =
                    meshOctreeCubeCoordinates::edgeNodes_[fEdges[feJ]];
                
                const label sp = sEdge[0];
                const label ep = eEdge[0];
                
                if( pl[sp] == pl[ep] )
                    continue;
                if( pl[sp] != pl[sEdge[1]] )
                    continue;
                if( pl[ep] != pl[eEdge[1]] )
                    continue;
                
                # ifdef DEBUGTets
                Info << "Octree node " << pointI << " has leaves";
                for(label plI=0;plI<8;++plI)
                    Info << ' ' << pointLeaves(pointI, plI);
                Info << endl;
                Info << "Searching face " << fI << endl;
                Info << "Searching face edge " << feI << endl;
                # endif
                
                //- allocate face centre labels
                if( faceCentreLabel.sizeOfRow(pl[sp]) == 0 )
                {
                    faceCentreLabel.setRowSize(pl[sp], 6);
                    forAllRow(faceCentreLabel, pl[sp], k)
                        faceCentreLabel(pl[sp], k) = -1;
                }
                if( faceCentreLabel.sizeOfRow(pl[ep]) == 0 )
                {
                    faceCentreLabel.setRowSize(pl[ep], 6);
                    forAllRow(faceCentreLabel, pl[ep], k)
                        faceCentreLabel(pl[ep], k) = -1;
                }
                //- create centre labels
                label fs, fe;
                
                fs = meshOctreeCubeCoordinates::edgeFaces_[fEdges[feJ]][0];
                if( fs == fI )
                    fs = meshOctreeCubeCoordinates::edgeFaces_[fEdges[feJ]][1];
                
                fe = meshOctreeCubeCoordinates::edgeFaces_[fEdges[feI]][0];
                if( fe == fI )
                    fe = meshOctreeCubeCoordinates::edgeFaces_[fEdges[feI]][1];
                
                # ifdef DEBUGTets
                Info << "Face for the cube at edge " << feI << " is "
                     << fs << endl;
                Info << "Face for the cube at edge " << feJ << " is "
                     << fe << endl;
                #endif
                
                //- create face centre point
                if( faceCentreLabel(pl[sp], fs) == -1 )
                {
                    faceCentreLabel(pl[sp], fs) = tetPoints_.size();
                    faceCentreLabel(pl[ep], fe) = tetPoints_.size();
                    point p(vector::zero);
                    const label* fn = meshOctreeCubeCoordinates::faceNodes_[fe];
                    for(label i=0;i<4;++i)
                        p += tetPoints_[nodeLabels(pl[ep], fn[i])];
                    p /= 4;
                    tetPoints_.append(p);
                }
                
                //- create tets connecting centroids
                checkAndAppendTet
                (
                    partTet
                    (
                        faceCentreLabel(pl[sp], fs),
                        subNodeLabels(pl[ep], 7-eEdge[0]),
                        subNodeLabels(pl[ep], 7-eEdge[1]),
                        cubeLabel[pl[ep]]
                    )
                );
                checkAndAppendTet
                (
                    partTet
                    (
                        faceCentreLabel(pl[sp], fs),
                        subNodeLabels(pl[sp], 7-sEdge[1]),
                        subNodeLabels(pl[sp], 7-sEdge[0]),
                        cubeLabel[pl[sp]]
                    )
                );
                
                FixedList<label, 4> subNodes;
                subNodes[0] = subNodeLabels(pl[sp], 7-sEdge[1]);
                subNodes[1] = subNodeLabels(pl[sp], 7-sEdge[0]);
                subNodes[2] = subNodeLabels(pl[ep], 7-eEdge[0]);
                subNodes[3] = subNodeLabels(pl[ep], 7-eEdge[1]);
                
                # ifdef DEBUGTets
                Info << "Sub nodes are " << subNodes << endl;
                # endif
                
                forAll(subNodes, nodeI)
                {
                    checkAndAppendTet
                    (
                        partTet
                        (
                            pointI,
                            subNodes[nodeI],
                            subNodes[(nodeI+1)%4],
                            faceCentreLabel(pl[sEdge[0]], fs)
                        )
                    );
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
