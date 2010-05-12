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
//#define debugCheck

namespace Foam
{

void binaryTree::insert
( 
    const chemPoint& x, 
    binaryNode *t
) 
{

#   ifdef debugCheck       
    Info << "binaryTree::insert::Checking the node" << endl;

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

    Info << "binaryTree::insert::End checking the node" << endl;
#   endif


//  Checking if the list is full


            
    //  First possibility. empty node. insert the element on the left side.
    /*
                O               O
               / \      ->     / \ 
              n   n           E   n
    */
    if
    (
        t->elementLeft() == NULL &&
        t->elementRight() == NULL &&
        t->left() == NULL &&
        t->right() == NULL
    )
    {
#       ifdef debugCheck       
        Info << "First possibility" << endl;
        Info << "element on the left" << endl;
#       endif
        insertElementLeft(x, t);
    }

    //  Second possibility. no element on the right side. insert the element on the right side.
    /*
                O               O
               / \      ->     / \  
              E   n           E   E
    */

    else if
    (
        t->elementLeft() != NULL &&
        t->elementRight() == NULL &&
        t->left() == NULL &&
        t->right() == NULL
    )
    {
#       ifdef debugCheck       
        Info << "Second possibility" << endl;
        Info << "element on the right" << endl;
#       endif
        insertElementRight(x, t);
    }

    //  Third possibility. Elements in both sides. Insert a new node.
    /*
                O               O                       O
               / \             / \                     / \      
              E   E     ->    E   O         OR        O   E 
                                 / \                 / \
                                E   E               E   E
    */

    else if
    ( 
        t->elementLeft() != NULL &&
        t->elementRight() != NULL &&
        t->left() == NULL &&
        t->right() == NULL
    )
    {
#       ifdef debugCheck       
        Info << "Third possibility" << endl;
        Info << "Full node, insert new node"  << endl;
#       endif
        insertNewNode(x, t);   
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
#       ifdef debugCheck       
        Info << "Fourth possibility" << endl;
#       endif
        // Check if the point is on the right
               
        if(onRight(x, t))
        {   
#           ifdef debugCheck       
            Info << "newNode on the right" << endl;
#           endif    
            insertNewNodeRight(x, t);   
        }
        else
        {
#           ifdef debugCheck       
            Info << "travel to the next node on the left" << endl;
#           endif
//            insert(x, t->left());
            insert(x, t->left_);
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
#       ifdef debugCheck       
        Info << "Fifth possibility" << endl;
#       endif
        if(onLeft(x, t))
        {
#           ifdef debugCheck       
            Info << "newNode on the left" << endl;
#           endif
            insertNewNodeLeft(x, t);   
        }
        else
        {
#           ifdef debugCheck       
            Info << "travel to the next node on the right" << endl;
#           endif
            insert(x, t->right_);
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
        
#       ifdef debugCheck       
        Info << "Sixth possibility" << endl;
        Info << "travel" << endl;
#       endif        
        if(onRight(x, t))
        {    
            insert(x, t->right());   
        }
        else
        {
            insert(x, t->left());   
        }
    }
    else
    {

        FatalError << "binaryTree::insert, the node has wrong pointers" << abort(FatalError);
#       ifdef debugCheck       

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

#       endif        

   }

#       ifdef debugCheck       
        Info << "Exiting" << endl;
#       endif        

        
}

void binaryTree::insertElementLeft
( 
    const chemPoint& x, 
    binaryNode *t
) 
{
    
    label n = chemPointList_.size();

    chemPointList_.setSize(n+1);
    chemPointList_.append(new chemPoint(x,new binaryNode()));  

    chemPoint *p = chemPointList_[n];

    chemPointList_[n]->node_ = t;
    chemPointList_[n]->node()->elementLeft_ = p;

    scalarField& vLeft = chemPointList_[n]->node()->vLeft();
    vLeft.setSize(p->v0().size());
    vLeft = scaleComposition(p->v0());
    
    
}

void binaryTree::insertElementRight
( 
    const chemPoint & x, 
    binaryNode *t
) 
{

    label n = chemPointList_.size();

    chemPointList_.setSize(n+1);
    chemPointList_.append(new chemPoint(x,new binaryNode()));  

    chemPoint *p = chemPointList_[n];

    chemPointList_[n]->node_ = t;
    chemPointList_[n]->node()->elementRight_ = p;

    scalarField& vRight = chemPointList_[n]->node()->vRight();
    vRight.setSize(p->v0().size());
    vRight = scaleComposition(p->v0());

#   include "calculatevTa.H"    

/*

    chemPointList_[n].node_.vT() = chemPointList_[n].node_.elementRight_->v0() - chemPointList_[n].node_.elementLeft_->v0();
    scalarField fiMedium = 0.5*(chemPointList_[n].node_.elementLeft_->v0() + chemPointList_[n].node_.elementRight_->v0());
    chemPointList_[n].node_.a() = sum(fiMedium*vT());  

*/

}

void binaryTree::insertNewNode
( 
    const chemPoint & x, 
    binaryNode *t
) 
{

    if(onLeft(x, t))
    {
        // element on the left side
        insertNewNodeLeft(x, t);  
    }
    else
    {
        // element on the right side 
        insertNewNodeRight(x, t);   
    }
    
}

void binaryTree::insertNewNodeLeft
(
    const chemPoint& x, 
    binaryNode *t
)
{

    if(isListFull())
    {

//        Info << "list full left" << endl;
        
        chemPoint &p = findLessUsed(nodeList_);
/*        
        if(p.node() == NULL)
        {
            Info << "p.node() = null" << endl;
        }
        if(p.node()->elementRight_ == NULL)
        {
            Info << "p.node().elementRight = null" << endl;
        }
        if(p.node()->elementLeft_ == NULL)
        {
            Info << "p.node().elementLeft = null" << endl;
        }
*/        
        Info << "LIST FULL LEFT" << endl;
        clean(p, p.node(), x, t);
        add(x, root_, p, p.node());
    }
    else
    {


#       ifdef debugCheck       
        Info << "Insert node on the left"<< endl;
#       endif
       
        label m = nodeList_.size();
        nodeList_.setSize(m+1);
        nodeList_.append(new binaryNode(NULL, NULL, NULL, NULL, NULL));
        binaryNode* node = nodeList_[m];
        
        label n = chemPointList_.size();       
        chemPointList_.setSize(n+1);
        chemPointList_.append(new chemPoint(x,node));
        chemPoint* cPoint = chemPointList_[n];


//        binaryNode* node = &nodeList_[m];
        t->left_ = node;
        
        node->down_ = t;
        node->elementLeft_ = t->elementLeft_;

        node->elementLeft_->node_ = node;

        t->elementLeft_ = NULL;        
        
//        node->elementRight_ = &chemPointList_[n];  
        node->elementRight_ = cPoint;  

//        node->elementLeft_->node_ = &nodeList_[m];
//        node->elementRight_->node_ = &nodeList_[m];
        node->elementLeft_-> node_ = node;
        node->elementRight_-> node_ = node;
        
/*
        node->vT() = node->elementRight_->v0() - node->elementLeft_->v0();

        scalarField fiMedium = 0.5*(node->elementLeft_->v0() + node->elementRight_->v0());
        node->a() = sum(fiMedium*node->vT());  
*/

        scalarField vsR = scaleComposition(node->elementRight_->v0());
        scalarField vsL = scaleComposition(node->elementLeft_->v0());
        
        node->vT().setSize(vsR.size());
        node->vT() = vsR - vsL;
        scalarField fiMedium = 0.5 * (vsR + vsL);
        node->a() = sum(fiMedium*node->vT());
        
        node->left_ = NULL;
        node->right_ = NULL;
        
        scalarField& vLeft = node->vLeft();
        vLeft.setSize(node->elementLeft_->v0().size());
        vLeft = scaleComposition(node->elementLeft_->v0());
                
        scalarField& vRight = node->vRight();
        vRight.setSize(node->elementRight_->v0().size());
        vRight = scaleComposition(node->elementRight_->v0());


/*
        Info << "node->vT()" << node->vT() << endl;

        Info << "node->a()" << node->a() << endl;
*/        
#       ifdef debugCheck       
        Info << "end insert node on the left"<< endl;
#       endif

    }

}

void binaryTree::insertNewNodeRight
(
    const chemPoint & x, 
    binaryNode *t
)
{

    if(isListFull())
    {

//        Info << "list full right" << endl;

        chemPoint &p = findLessUsed(nodeList_);
/*        
        if(p.node() == NULL)
        {
            Info << "p.node() = null" << endl;
        }
        if(p.node()->elementRight_ == NULL)
        {
            Info << "p.node().elementRight = null" << endl;
        }
        if(p.node()->elementLeft_ == NULL)
        {
            Info << "p.node().elementLeft = null" << endl;
        }
*/        
        Info << "LIST FULL RIGHT" << endl;
        clean(p, p.node(), x, t);
        add(x, root_, p, p.node());
    }
    else
    {

#       ifdef debugCheck       
        Info << "Insert node on the right"<< endl;
#       endif
       
        label m = nodeList_.size();
        nodeList_.setSize(m+1);
        nodeList_.append(new binaryNode(NULL, NULL, NULL, NULL, NULL));
        binaryNode *node = nodeList_[m];
        
        label n = chemPointList_.size();       
        chemPointList_.setSize(n+1);
        chemPointList_.append(new chemPoint(x,node));
        chemPoint* cPoint = chemPointList_[n];

//        binaryNode *node = &nodeList_[m];
        t->right_ = node;
        
        node->down_ = t;
        node->elementLeft_ = t->elementRight_;        
        t->elementRight_ = NULL;        
        
//        node->elementRight_ = &chemPointList_[n];        
        node->elementRight_ = cPoint;        

//        node->elementLeft_->node_ = &nodeList_[m];
//        node->elementRight_->node_ = &nodeList_[m];
        node->elementLeft_->node_ = node;
        node->elementRight_->node_ = node;


/*
        node->vT() = node->elementRight_->v0() - node->elementLeft_->v0();

        scalarField fiMedium = 0.5*(node->elementLeft_->v0() + node->elementRight_->v0());
        node->a() = sum(fiMedium*node->vT());  
*/

        scalarField vsR = scaleComposition(node->elementRight_->v0());
        scalarField vsL = scaleComposition(node->elementLeft_->v0());
                
        node->vT() = vsR - vsL;
        node->vT().setSize(vsR.size());
        node->vT() = vsR - vsL;
        
        
        scalarField fiMedium = 0.5 * (vsR + vsL);
        node->a() = sum(fiMedium*node->vT());
        
        node->left_ = NULL;
        node->right_ = NULL;

        node->vLeft() = scaleComposition(node->elementLeft_->v0());
        node->vRight() = scaleComposition(node->elementRight_->v0());
    
        scalarField& vLeft = node->vLeft();
        vLeft.setSize(node->elementLeft_->v0().size());
        vLeft = scaleComposition(node->elementLeft_->v0());
                
        scalarField& vRight = node->vRight();
        vRight.setSize(node->elementRight_->v0().size());
        vRight = scaleComposition(node->elementRight_->v0());


/*
        Info << "node->vT()" << node->vT() << endl;

        Info << "node->a()" << node->a() << endl;
*/        

#       ifdef debugCheck       
        Info << "end insert node on the right"<< endl;
#       endif
    
    }

}


}
