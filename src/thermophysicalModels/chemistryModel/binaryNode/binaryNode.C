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

#include "chemPoint.H"
#include "scalarField.H"
#include "binaryNode.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

binaryNode::binaryNode
(
)
:
    elementLeft_(NULL),
    elementRight_(NULL),
    left_(NULL), 
    right_(NULL),
    down_(NULL),
    vLeft_(0,0),
    vRight_(0,0)
{}


binaryNode::binaryNode
(
    chemPoint* elementLeft,
    chemPoint* elementRight,
    binaryNode* left,
    binaryNode* right,
    binaryNode* down
)
:
    elementLeft_(elementLeft),
    elementRight_(elementRight),
    left_(left), 
    right_(right),
    down_(down),
    vLeft_(0,0),
    vRight_(0,0)
{}


binaryNode::binaryNode
(
    binaryNode *bn
)
:
    elementLeft_(bn->elementLeft()),
    elementRight_(bn->elementRight()),
    left_(bn->left()), 
    right_(bn->right()),
    down_(bn->down()),
    vLeft_(0,0),
    vRight_(0,0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void binaryNode::setFree()
{
    vLeft_.setSize(0);
    vRight_.setSize(0);
    vT_ = 0.0;
    a_ = 0.0;
    

    elementLeft_ = NULL;
    elementRight_ = NULL;
    
    
    left_ = NULL ;
    right_ = NULL;
    down_ = NULL;
}


void binaryNode::updateNode()
{
    vT_ = vRight_ - vLeft_;
    scalarField fiMedium = 0.5*(vRight_ + vLeft_);
    a_ = sum(fiMedium*vT_);
}


void binaryNode::clearData()
{

//        Info << "clear data" << endl;

/*        
        if(vT_.size()>0)
        {
            Info << "vT" << endl;
            vT_.clear();
        }
        if(vLeft_.size()>0)
        {
            Info << "vleft" << endl;
            vLeft_.clear();
        }
        if(vRight_.size()>0)
        {
            Info << "vright" << endl;
            vRight_.clear();
        }
        
        Info << "clear()" << endl;

*/        

/*
        if(left_)
        {
            Info << "left" << endl;
//            left_ = NULL ;
        }
        if(right_)
        {
            Info << "right" << endl;
//            right_ = NULL;
        }
        if(down_)
        {
            Info << "down" << endl;
//            down_ = NULL;
        }

*/

        if(elementLeft_)
        {
//            Info << "elementLeft" << endl;
            elementLeft_->clearData();
//            if(elementLeft_->node_)
            {
//                elementLeft_->node_=NULL;
            }
            deleteDemandDrivenData(elementLeft_);
        }
        
        if(elementRight_)
        {
//            Info << "elementRight" << endl;
            elementRight_->clearData();
//            if(elementRight_->node_)
            {
//                elementRight_->node_=NULL;
            }
            deleteDemandDrivenData(elementRight_);
        }
//        Info << "deleted" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
