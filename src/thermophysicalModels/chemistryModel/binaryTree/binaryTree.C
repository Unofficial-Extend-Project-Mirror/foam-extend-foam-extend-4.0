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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "binaryTree.H"
#include "binaryNode.H"
#include "chemistryOnlineLibrary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

binaryTree::binaryTree
(
    const dictionary treeDict,
    chemistryOnlineLibrary& onlineLibrary
) 
:
    start_(NULL, NULL, NULL, NULL, NULL),
    maxElements_(readLabel(treeDict.lookup("maxElements"))),
    tolerance_(readScalar(treeDict.lookup("tolerance"))),
    onlineLibrary_(onlineLibrary)
    {
#       ifdef debugCheck       
        Info << "binaryTree::binaryTree -> root assignment" << endl;
#       endif

        root_ = &start_;

#       ifdef debugCheck       
        Info << "binaryTree::binaryTree -> end root assignment" << endl;
#       endif
    }


scalarField binaryTree::scaleComposition(const scalarField& x) const
{
            
    scalarField xScaled = x;
    scalarField scaleF = onlineLibrary_.scaleFactor();
    
    for(label i=0; i < xScaled.size(); i++)
    {
        if(scaleF[i] < VSMALL)
        {
            xScaled[i] = 0;
        }
        else
        {
            xScaled[i] /= scaleF[i];
        }
    }    

    return xScaled;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


