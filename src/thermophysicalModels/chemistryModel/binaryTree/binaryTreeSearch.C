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
//#define debugSearch

namespace Foam
{


chemPoint* binaryTree::findClosest
(
    const scalarField& x, 
    binaryNode *t
)
{
    

//  Exceptions to deal with:
//
//  1. Empty Node
//  2. Node with one element (on the left...)

#       ifdef debugSearch
        Info << "binaryTree::findClosest" << endl;
#       endif
    
    if
    (
        t->elementLeft() == NULL &&
        t->elementRight() == NULL &&
        t->left() == NULL &&
        t->right() == NULL
    )
    {
#       ifdef debugSearch
        Info << "Exception 1" << endl;
#       endif
        return NULL;
    }
    
    if
    (
        t->elementLeft() != NULL &&
        t->elementRight() == NULL &&
        t->left() == NULL &&
        t->right() == NULL   
    )
    {

#       ifdef debugSearch
        Info << "Exception 2" << endl;
#       endif
/*
    
        if(t->elementLeft_->difference(x) < tolerance_)
        {
#           ifdef debugSearch
            Info << "FOUND on the ELEMENT LEFT!!!" << endl;
            Info << "Difference = " << t->elementRight_->difference(x) << endl;
#           endif
            return t->elementLeft_;
        }
        else
        {
            return NULL;
        }   
*/
        return t->elementLeft_;

    }
    
    
       
    if(onRight(x,t))
    {
    
        if(t->right_ != NULL)  // 8-)
        {

#           ifdef debugSearch
            Info << "Travel on right!!!" << endl;
#           endif
        
            return findClosest(x,t->right_);
            
        }
        else if(t->elementRight_ != NULL)  
        {
        

/*
            if(t->elementRight_->difference(x) < tolerance_)
            {
#               ifdef debugSearch
                Info << "FOUND on the ELEMENT RIGHT!!!" << endl;
                Info << "Difference = " << t->elementRight_->difference(x) << endl;
#               endif
                return t->elementRight_;
            }
            else
            {
                return NULL;
            }
*/
            return t->elementRight_;
        
        }
        else    //the node is empty
        {

            return NULL;

        }
        
    }
    else
    {
    
        if(t->left_ != NULL)  // 8-)
        {

#           ifdef debugSearch
            Info << "Travel on left!!!" << endl;
#           endif
        
            return findClosest(x,t->left_);

        }
        else if(t->elementLeft_ != NULL)
        {

/*
            if(t->elementLeft_->difference(x) < tolerance_)
            {
#               ifdef debugSearch
                Info << "FOUND on the ELEMENT LEFT!!!" << endl;
                Info << "Difference = " << t->elementLeft_->difference(x) << endl;
#               endif
                return t->elementLeft_;
            }
            else
            {
                return NULL;
            }
*/

            return t->elementLeft_;

        }
        else    //the node is empty
        {
            return NULL;
        }
        
    }        
    
    return NULL;    

}



}
