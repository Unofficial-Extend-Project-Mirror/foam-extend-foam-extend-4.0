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
//#define debugAdd
//#define debugPreAdd

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void binaryTree::add
( 
    const chemPoint& x, 
    binaryNode *t,
    chemPoint& oldPoint,
    binaryNode *oldNode
) 
{

        #   ifdef debugAdd       
        Info << "binaryTree::add::Checking the node" << endl;

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

        Info << "binaryTree::add::End checking the node" << endl;
        #   endif
    
   

    //  Third possibility. Elements in both sides. Insert a new node.
    /*
                O               O                       O
               / \             / \                     / \      
              E   E     ->    E   O         OR        O   E 
                                 / \                 / \
                                E   E               E   E
    */

    if
    ( 
        t->elementLeft() != NULL &&
        t->elementRight() != NULL &&
        t->left() == NULL &&
        t->right() == NULL
    )
    {
        addNewNode(x, t, oldPoint, oldNode);   
    }

    //  Fourth possibility. Elements on the left, node on the right. What can I do?.
    /*
                O                O               
               / \              / \              
              O   E      ->    O   O        or to travel...    
                                  / \            
                                 E   E           
    */

    else if
    (
        t->elementLeft() == NULL &&
        t->elementRight() != NULL &&
        t->left() != NULL &&
        t->right() == NULL
    )
    {
        // Check if the point is on the right
               
        if(onRight(x, t))
        {   
            addNewNodeRight(x, t, oldPoint, oldNode);   
        }
        else
        {
            add(x, t->left_, oldPoint, oldNode);
        }
        
    }
    //  Fifth possibility. Elements on the right, node on the left. Insert a new node on the left side.
    /*
                O                O    
               / \              / \   
              E   O      ->    O   O        or  to travel...    
                              / \     
                             E   E    
    */

    else if
    (
        t->elementLeft() != NULL &&
        t->elementRight() == NULL &&
        t->left() == NULL &&
        t->right() != NULL
    )
    {
        if(onLeft(x, t))
        {
            addNewNodeLeft(x, t, oldPoint, oldNode);   
        }
        else
        {
            add(x, t->right_, oldPoint, oldNode);
        }
    }

    //  Sixth possibility. Nodes on both sides. Go to the next node.
    /*
                O
               / \
              O   O         
    */

    else if
    (
        t->elementLeft() == NULL &&
        t->elementRight() == NULL &&
        t->left() != NULL &&
        t->right() != NULL
    )
    {
        
        if(onRight(x, t))
        {    
            add(x, t->right(), oldPoint, oldNode);   
        }
        else
        {
            add(x, t->left(), oldPoint, oldNode);   
        }
    }
        
}



void binaryTree::addNewNode
( 
    const chemPoint & x, 
    binaryNode *t,
    chemPoint& oldPoint,
    binaryNode *oldNode
) 
{

#   ifdef debugAdd       
    Info << "binaryTree::addNewNode" << endl;
#   endif    
    
    if(onLeft(x, t))
    {
        // element on the left side
#   ifdef debugAdd       
    Info << "binaryTree::addNewNode::addNewNodeLeft" << endl;
#   endif    
        addNewNodeLeft(x, t, oldPoint, oldNode);  
    }
    else
    {
        // element on the right side 
#   ifdef debugAdd       
    Info << "binaryTree::addNewNode::addNewNodeRight" << endl;
#   endif    
        addNewNodeRight(x, t, oldPoint, oldNode);   
    }
    
}

void binaryTree::addNewNodeLeft
(
    const chemPoint& x, 
    binaryNode *t,
    chemPoint& oldPoint,
    binaryNode *oldNode
)
{


        #   ifdef debugPreAdd       
        Info << "binaryTree::addNewNodeLeft::Checking the parente node BEFORE" << endl;

        if
        (
            t->elementLeft() == NULL 
        )
        {
            Info << "t->elementLeft() == NULL" << endl;   
        }
        if
        (
            t->elementRight() == NULL
        )
        {
            Info << "t->elementRight() == NULL" << endl;   
        }
        if
        (
            t->left() == NULL
        )
        {
            Info << "t->left() == NULL" << endl;       
        }
        if
        (
            t->right() == NULL
        )
        {
            Info << "t->right() == NULL" << endl;       
        }

        Info << "binaryTree::addNewNodeLeft::End checking the parent node BEFORE" << endl;
        #   endif

        oldPoint.setFree();
        oldNode->setFree();
 
//        oldPoint = chemPoint(x,new binaryNode());
        oldPoint = chemPoint(x,oldNode);

        t->left_ = oldNode;

        oldNode->down_ = t;
        oldNode->left_ = NULL;
        oldNode->right_ = NULL;
       
        oldNode->elementLeft_ = t->elementLeft_;
        oldNode->elementLeft_->node_ = oldNode;
            
        t->elementLeft_ = NULL;  

        oldNode->elementRight_ = &oldPoint;  
                
        scalarField vsR = scaleComposition(oldNode->elementRight_->v0());
        scalarField vsL = scaleComposition(oldNode->elementLeft_->v0());

//        oldNode->vT() = oldNode->elementRight_->v0() - oldNode->elementLeft_->v0();
        oldNode->vT() = vsR - vsL;

//        scalarField fiMedium = 0.5*(oldNode->elementLeft_->v0() + oldNode->elementRight_->v0());
        scalarField fiMedium = 0.5*(vsR + vsL);
        oldNode->a() = sum(fiMedium*oldNode->vT());  

        oldPoint.node_= oldNode;
        
        scalarField& vLeft = oldNode->vLeft();
        vLeft.setSize(oldNode->elementLeft_->v0().size());
        vLeft = scaleComposition(oldNode->elementLeft_->v0());

        scalarField& vRight = oldNode->vRight();
        vRight.setSize(oldNode->elementRight_->v0().size());
        vRight = scaleComposition(oldNode->elementRight_->v0());

        #   ifdef debugAdd       
        Info << "binaryTree::addNewNodeLeft::Checking the node" << endl;

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

        Info << "binaryTree::addNewNodeLeft::End checking the node" << endl;
        #   endif

        #   ifdef debugAdd       
        Info << "binaryTree::addNewNodeLeft::Checking the PARENT node" << endl;

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

        Info << "binaryTree::addNewNodeLeft::End checking the PARENT node" << endl;
        #   endif

        

}

void binaryTree::addNewNodeRight
(
    const chemPoint & x, 
    binaryNode *t,
    chemPoint& oldPoint,
    binaryNode *oldNode
)
{
#   ifdef debugPreAdd

    Info << "binaryTree::addNewNodeRight::Checking the parente node BEFORE"
        << endl;

    if
    (
        t->elementLeft() == NULL 
    )
    {
        Info << "t->elementLeft() == NULL" << endl;   
    }
    if
    (
        t->elementRight() == NULL
    )
    {
        Info << "t->elementRight() == NULL" << endl;   
    }
    if
    (
        t->left() == NULL
    )
    {
        Info << "t->left() == NULL" << endl;       
    }
    if
    (
        t->right() == NULL
    )
    {
        Info << "t->right() == NULL" << endl;       
    }

    Info<< "binaryTree::addNewNodeRight::End checking the parent node BEFORE"
        << endl;
#   endif

    oldPoint.setFree();
    oldNode->setFree();        

//     oldPoint = chemPoint(x,new binaryNode());
    oldPoint = chemPoint(x,oldNode);

    t->right_ = oldNode;

    oldNode->down_ = t;
    oldNode->left_ = NULL;
    oldNode->right_ = NULL;

    oldPoint.node_= oldNode;

//    Info << t->elementRight_->v0() << endl;

    t->elementRight_->node_ = oldNode;

    oldNode->elementLeft_ = t->elementRight_;

//    oldNode->elementLeft_ = t->elementRight_;        
//    oldNode->elementLeft_->node_ = oldNode;        
        
    t->elementRight_ = NULL;        

    oldNode->elementRight_ = &oldPoint;        
        
/*
    oldNode->vT() = oldNode->elementRight_->v0() - oldNode->elementLeft_->v0();

    scalarField fiMedium = 0.5*(oldNode->elementLeft_->v0() + oldNode->elementRight_->v0());

    oldNode->a() = sum(fiMedium*oldNode->vT());  
*/

    scalarField vsR = scaleComposition(oldNode->elementRight_->v0());
    scalarField vsL = scaleComposition(oldNode->elementLeft_->v0());

    oldNode->vT() = vsR - vsL;

    scalarField fiMedium = 0.5*(vsR + vsL);
    oldNode->a() = sum(fiMedium*oldNode->vT());  


    oldNode->left_ = NULL;
    oldNode->right_ = NULL;

    scalarField& vLeft = oldNode->vLeft();
    vLeft.setSize(oldNode->elementLeft_->v0().size());
    vLeft = scaleComposition(oldNode->elementLeft_->v0());

    scalarField& vRight = oldNode->vRight();
    vRight.setSize(oldNode->elementRight_->v0().size());
    vRight = scaleComposition(oldNode->elementRight_->v0());


//    oldPoint.node_= oldNode;


/*
//      NEW!!!

    t->right_ = oldNode;

    oldNode->down_ = t;
    oldNode->elementLeft_ = t->elementRight_;        
    t->elementRight_ = NULL;        

    oldNode->elementRight_ = &oldPoint;        

    oldNode->elementLeft_->node_ = oldNode;
    oldNode->elementRight_->node_ = oldNode;

    oldNode->vT() = oldNode->elementRight_->v0() - oldNode->elementLeft_->v0();

    scalarField fiMedium =
    0.5*(oldNode->elementLeft_->v0() + oldNode->elementRight_->v0());

    oldNode->a() = sum(fiMedium*oldNode->vT());  

    oldNode->left_ = NULL;
    oldNode->right_ = NULL;

*/

#   ifdef debugAdd       
    Info << "binaryTree::addNewNodeRight::Checking the node" << endl;

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

    Info << "binaryTree::addNewNodeRight::End checking the node" << endl;
#   endif
   

#   ifdef debugAdd       
    Info << "binaryTree::addNewNodeRight::Checking the PARENT node" << endl;

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

    Info << "binaryTree::addNewNodeRight::End checking the PARENT node" << endl;
#   endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
