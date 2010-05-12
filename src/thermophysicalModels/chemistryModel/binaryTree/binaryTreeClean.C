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
#include "demandDrivenData.H"
//#define debugFindLessUsed
//#define debugReplace

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool binaryTree::isListFull()
{
    label n = chemPointList_.size();

#   ifdef debugReplace       
    Info << "n " << n << endl;
    Info << "maxElements " << maxElements_ << endl;
#   endif       
    
    if(n <= maxElements_)
    {
        return false;
    }
    else
    {   
        return true;
    }
      
}


bool binaryTree::ListFull()
{
    return isListFull();
}


chemPoint& binaryTree::findLessUsed
(
    DynamicList<chemPoint*> pointList
)
{
    
#   ifdef debugFindLessUsed
    Info << "binaryTree::findLessUsed(PtrList<chemPoint>& pointList)" << endl;
#   endif
    
    label nLessUsed = pointList[0]->nUsed();
    
    Info << "nLessUsed" << endl;
    
    label leafLessUsed = 0;
    forAll(pointList, pointi)
    {
        if(pointList[pointi]->nUsed() < nLessUsed)
        {
            nLessUsed = pointList[pointi]->nUsed();
            Info << "nLessUsed" << endl;
            leafLessUsed = pointi;
            Info << "nLessUsed" << endl;
        }
    }
    
    return *pointList[leafLessUsed];

}

chemPoint& binaryTree::findLessUsed
(
    DynamicList<binaryNode*> nodeList
)
{
    
   
//    label nLessUsed = (label)VGREAT;
//    label nodeLessUsed = (label)VGREAT;
    
    label nLessUsed = INT_MAX;
    label nodeLessUsed = INT_MAX;
    
    word side = "left";

    forAll(nodeList, nodei)
    {
        
        if(nodeList[nodei]->elementLeft_ != NULL)
        {
        
            if(nodeList[nodei]->elementLeft_->nUsed() <= nLessUsed)
            {
                nLessUsed = nodeList[nodei]->elementLeft_->nUsed();
                nodeLessUsed = nodei;
                side = "left";
#               ifdef debugFindLessUsed
                Info << nodeList[nodei]->elementLeft_->v0() << endl;
                Info << nodeList[nodei]->elementLeft_->nUsed() << endl;
#               endif
            }
        }
        
        if(nodeList[nodei]->elementRight_ != NULL)
        {
            if(nodeList[nodei]->elementRight_->nUsed() <= nLessUsed)
            {
                nLessUsed = nodeList[nodei]->elementRight_->nUsed();
                nodeLessUsed = nodei;
                side = "right";
#               ifdef debugFindLessUsed
                Info << nodeList[nodei]->elementRight_->v0() << endl;
                Info << nodeList[nodei]->elementRight_->nUsed() << endl;
#               endif
            }
        }
    }
    
    
    if(side == "right")
    {
#       ifdef debugFindLessUsed
        Info  << "value RIGHT "<< nodeList[nodeLessUsed]->elementRight_->v0() << endl;
#       endif
        return *nodeList[nodeLessUsed]->elementRight_;
    }
    else
    {
#       ifdef debugFindLessUsed
        Info << "value LEFT " << nodeList[nodeLessUsed]->elementLeft_->v0() << endl;
#       endif
        return *nodeList[nodeLessUsed]->elementLeft_;
    }
}

chemPoint& binaryTree::findLessUsed
(
    const chemPoint& p,
    DynamicList<binaryNode*> nodeList
)
{
    
   
//    label nLessUsed = (label)VGREAT;
//    label nodeLessUsed = (label)VGREAT;
    
    label nLessUsed = INT_MAX;
    label nodeLessUsed = INT_MAX;

//  p is the point to insert
//    scalar norm = 0.0;
    
    word side = "left";

    forAll(nodeList, nodei)
    {
        
        if(nodeList[nodei]->elementLeft_ != NULL)
        {
        
            if(nodeList[nodei]->elementLeft_->nUsed() <= nLessUsed)
            {
                nLessUsed = nodeList[nodei]->elementLeft_->nUsed();
                nodeLessUsed = nodei;
                side = "left";
#               ifdef debugFindLessUsed
                Info << nodeList[nodei]->elementLeft_->v0() << endl;
                Info << nodeList[nodei]->elementLeft_->nUsed() << endl;
#               endif
            }
        }
        
        if(nodeList[nodei]->elementRight_ != NULL)
        {
            if(nodeList[nodei]->elementRight_->nUsed() <= nLessUsed)
            {
                nLessUsed = nodeList[nodei]->elementRight_->nUsed();
                nodeLessUsed = nodei;
                side = "right";
#               ifdef debugFindLessUsed
                Info << nodeList[nodei]->elementRight_->v0() << endl;
                Info << nodeList[nodei]->elementRight_->nUsed() << endl;
#               endif
            }
        }
    }
    
    
    if(side == "right")
    {
#       ifdef debugFindLessUsed
        Info  << "value RIGHT "<< nodeList[nodeLessUsed]->elementRight_->v0() << endl;
#       endif
        return *nodeList[nodeLessUsed]->elementRight_;
    }
    else
    {
#       ifdef debugFindLessUsed
        Info << "value LEFT " << nodeList[nodeLessUsed]->elementLeft_->v0() << endl;
#       endif
        return *nodeList[nodeLessUsed]->elementLeft_;
    }
}

void binaryTree::clean
(
    chemPoint& oldPoint, 
    binaryNode *oldNode, 
    const chemPoint& newPoint,
    binaryNode *parentNode
)
{
    
// 8-) 

#   ifdef debugReplace
    Info << "binaryTree::clean" << endl;
#   endif

   
#   ifdef debugReplace       
    Info << "binaryTree::clean::Checking the node" << endl;

    if
    (
        oldNode->elementLeft() == NULL 
    )
    {
        Info << "oldNode->elementLeft() == NULL" << endl;   
    }
    if
    (
        oldNode->elementRight() == NULL
    )
    {
        Info << "oldNode->elementRight() == NULL" << endl;   
    }
    if
    (
        oldNode->left() == NULL
    )
    {
        Info << "oldNode->left() == NULL" << endl;       
    }
    if
    (
        oldNode->right() == NULL
    )
    {
        Info << "oldNode->right() == NULL" << endl;       
    }

    Info << "binaryTree::clean::End checking the node" << endl;
#   endif

#   ifdef debugReplace       
    Info << "binaryTree::clean::Checking the PARENT node BEFORE" << endl;

    if
    (
        oldNode->down_->elementLeft() == NULL 
    )
    {
        Info << "oldNode->down_->elementLeft() == NULL" << endl;   
    }
    if
    (
        oldNode->down_->elementRight() == NULL
    )
    {
        Info << "oldNode->down_->elementRight() == NULL" << endl;   
    }
    if
    (
        oldNode->down_->left() == NULL
    )
    {
        Info << "oldNode->down_->left() == NULL" << endl;       
    }
    if
    (
        oldNode->down_->right() == NULL
    )
    {
        Info << "oldNode->down_->right() == NULL" << endl;       
    }

    Info << "binaryTree::clean::End checking the PARENT node BEFORE" << endl;
#   endif

//  Replace: connect the parent node with the following one, replacing the old node values with the new ones.

/*
        O
       / \          termination node - >       
      E   E
*/ 

    if(oldNode->elementRight() != NULL && oldNode->elementLeft() != NULL)
    {

#       ifdef debugReplace
        Info << "termination node" << endl;
#       endif

        //check if the element is on the right or on the left
        if(onRight(oldPoint,oldNode))
        {   
        
#           ifdef debugReplace
            Info << "node to replace on the right" << endl;
#           endif

            // The point to remove is on the right or on the left of the down node
            if(onRight(oldPoint,oldNode->down()))
            {
            
#               ifdef debugReplace
                Info << "node to replace on the right of the parent node" << endl;
#               endif

                //  The node to replace is on the right
                oldNode->down_->elementRight_ = oldNode->elementLeft_;
                oldNode->down_->elementRight_->node_ = oldNode->down_;
                oldNode->down_->right_ = NULL;

                //  Recalculating now vT and a for the oldNode->down



                scalarField& vRight = oldNode->down_->vRight();               
                vRight.setSize(oldNode->down_->elementRight_->v0().size());
                vRight = scaleComposition(oldNode->down_->elementRight_->v0());
                oldNode->down_->updateNode();

            }
            else
            {
            
#               ifdef debugReplace
                Info << "node to replace on the left of the parent node" << endl;
#               endif

                //  The node to replace is on the left
                oldNode->down_->elementLeft_ = oldNode->elementLeft_;           
                oldNode->down_->elementLeft_->node_ = oldNode->down_;
                oldNode->down_->left_ = NULL;

                //  Recalculating now vT and a for the oldNode->down

                scalarField& vLeft = oldNode->down_->vLeft();
                vLeft.setSize(oldNode->down_->elementLeft_->v0().size());
                vLeft = scaleComposition(oldNode->down_->elementLeft_->v0());
                oldNode->down_->updateNode();
                              
            }
            
        }   
        else
        {   
        
#           ifdef debugReplace
            Info << "node to replace on the left" << endl;
#           endif

            // The point to remove is on the left
            if(onRight(oldPoint,oldNode->down()))
            {
            
#               ifdef debugReplace
                Info << "node to replace on the right of the parent node" << endl;
#               endif

                //  The node to replace is on the right
                oldNode->down_->elementRight_ = oldNode->elementRight_;
                oldNode->down_->elementRight_->node_ = oldNode->down_;
                oldNode->down_->right_ = NULL;

                //  Recalculating now vT and a for the oldNode->down

                scalarField& vRight = oldNode->down_->vRight();               
                vRight.setSize(oldNode->down_->elementRight_->v0().size());
                vRight = scaleComposition(oldNode->down_->elementRight_->v0());
                oldNode->down_->updateNode();

            }
            else
            {
            
#               ifdef debugReplace
                Info << "node to replace on the left of the parent node" << endl;
#               endif

                //  The node to replace is on the right
                oldNode->down_->elementLeft_ = oldNode->elementRight_;           
                oldNode->down_->elementLeft_->node_ = oldNode->down_;
                oldNode->down_->left_ = NULL;

                //  Recalculating now vT and a for the oldNode->down

                scalarField& vLeft = oldNode->down_->vLeft();
                vLeft.setSize(oldNode->down_->elementLeft_->v0().size());
                vLeft = scaleComposition(oldNode->down_->elementLeft_->v0());

                oldNode->down_->updateNode();

            }               
        }
                  
    }
/*
        O                                      
       / \                                       O         
      E   O         intermediate node - >       / \         
         / \                                   .   .    
        .   .                                            
*/   
    else if(oldNode->right() != NULL && oldNode->elementLeft() != NULL)
    {
            
#       ifdef debugReplace
        Info << "Intermediate node: to replace on the left" << endl;
//        Info << "value" << oldPoint.v0() << endl;
#       endif

        oldNode->right_->down_ = oldNode->down_; 

        // connect the new parent node with the new node
        if(onRight(oldPoint, oldNode->down()))
        {
            oldNode->down_->right_ = oldNode->right_;
//            oldNode->right_->down_->right_ = oldNode->right_;
        }
        else
        {
            oldNode->down_->left_ = oldNode->right_;
        }       
        
    }
/*
        O
       / \                                      O 
      O   E        intermediate node - >       / \  
     / \                                      .   .   
    .   .         
*/   
    else if(oldNode->left()!= NULL  && oldNode->elementRight() != NULL)
    {

#       ifdef debugReplace
        Info << "Intermediate node: to replace on the right" << endl;
//        Info << "value" << oldPoint.v0() << endl;
#       endif

        oldNode->left_->down_ = oldNode->down_; 

        // connect the new parent node with the new node
        if(onRight(oldPoint, oldNode->down()))
        {
//            oldNode->left_->down_->right_ = oldNode->left_;
            oldNode->down_->right_ = oldNode->left_;
       }
        else
        {
//            oldNode->left_->down_->left_ = oldNode->left_;
            oldNode->down_->left_ = oldNode->left_;
        }
               
    }
    
    else
    {
//#       ifdef debugReplace
        FatalError << "binaryTree::clean, the node has wrong pointers" << abort(FatalError);
//#       endif        
    }
    
#   ifdef debugReplace
    Info << "binaryTree::clean" << endl;
#   endif

//  Replace: connect the parent node with the following one, replacing the old node values with the new ones.

/*
        O
       / \          termination node - >       
      E   E
*/ 
   
#   ifdef debugReplace       
    Info << "binaryTree::clean::Checking the node after replacement" << endl;

    if
    (
        oldNode->elementLeft() == NULL 
    )
    {
        Info << "oldNode->elementLeft() == NULL" << endl;   
    }
    if
    (
        oldNode->elementRight() == NULL
    )
    {
        Info << "oldNode->elementRight() == NULL" << endl;   
    }
    if
    (
        oldNode->left() == NULL
    )
    {
        Info << "oldNode->left() == NULL" << endl;       
    }
    if
    (
        oldNode->right() == NULL
    )
    {
        Info << "oldNode->right() == NULL" << endl;       
    }

    Info << "binaryTree::clean::End checking the node" << endl;
#   endif


//  End replacing nodes. Now replacing points. oldNode will become the "new" one    

#   ifdef debugReplace       
    Info << "binaryTree::clean::Checking the parent node" << endl;

    if
    (
        oldNode->down_->elementLeft() == NULL 
    )
    {
        Info << "oldNode->down_elementLeft() == NULL" << endl;   
    }
    if
    (
        oldNode->down_->elementRight() == NULL
    )
    {
        Info << "oldNode->down_->elementRight() == NULL" << endl;   
    }
    if
    (
        oldNode->down_->left() == NULL
    )
    {
        Info << "oldNode->down_->left() == NULL" << endl;       
    }
    if
    (
        oldNode->down_->right() == NULL
    )
    {
        Info << "oldNode->down_->right() == NULL" << endl;       
    }

    Info << "binaryTree::clean::End checking the parent node" << endl;
#   endif

}

void binaryTree::cleanAll()
{

//    deleteDemandDrivenData(root_->elementLeft_);
//    deleteDemandDrivenData(root_->elementRight_);
//    deleteDemandDrivenData(root_->right_);
//    deleteDemandDrivenData(root_->elementLeft_);

    root_->elementLeft_ = NULL;
    root_->elementRight_ = NULL;
    root_->right_ = NULL;
    root_->left_ = NULL;

    Info<< "Cleaning Chemistry Library: number of nodes "
        << nodeList_.size() << nl
        << "Cleaning Chemistry Library: number of points "
        << chemPointList_.size() << endl;

    forAll(nodeList_, nodeI)
    {
//        nodeList_[nodeI].clearData();
        deleteDemandDrivenData(nodeList_[nodeI]);
    }

    nodeList_.clear();

    forAll(chemPointList_, pointI)
    {
//        chemPoint cp(chemPointList_[pointI]);
//        chemPointList_[pointI].clearData();
//        chemPointList_[pointI] = cp;
        deleteDemandDrivenData(chemPointList_[pointI]) ;   
    }

    Info << "clearing chemPointList" << endl;           
    chemPointList_.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


}
