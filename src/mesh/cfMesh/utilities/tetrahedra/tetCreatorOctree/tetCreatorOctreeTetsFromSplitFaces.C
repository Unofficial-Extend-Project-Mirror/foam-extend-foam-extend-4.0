/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "tetCreatorOctree.H"
#include "demandDrivenData.H"

//#define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetCreatorOctree::createTetsFromSplitFaces()
{
    Info << "Creating tets from split faces" << endl;

    const labelList& cubeLabel = *cubeLabelPtr_;
    const VRWGraph& subNodeLabels = *subNodeLabelsPtr_;
    const FRWGraph<label, 8>& pLeaves = octreeCheck_.nodeLeaves();

    forAll(pLeaves, pointI)
    {
        for(label i=0;i<6;++i)
        {
            const label* fNodes =
                meshOctreeCubeCoordinates::faceNodes_[i];

            const label cLabel = pLeaves(pointI, fNodes[0]);

            if( cLabel < 0 )
                continue;
            if( cubeLabel[cLabel] < 0 )
                continue;

            if(
                (cLabel == pLeaves(pointI, fNodes[1])) &&
                (cLabel == pLeaves(pointI, fNodes[2])) &&
                (cLabel == pLeaves(pointI, fNodes[3]))
            )
            {
                //- create 4 tets
                for(label j=0;j<4;++j)
                {
                    checkAndAppendTet
                    (
                        partTet
                        (
                            pointI,
                            subNodeLabels(cLabel, 7-fNodes[j]),
                            subNodeLabels(cLabel, 7-fNodes[(j+1)%4]),
                            cubeLabel[cLabel]
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
