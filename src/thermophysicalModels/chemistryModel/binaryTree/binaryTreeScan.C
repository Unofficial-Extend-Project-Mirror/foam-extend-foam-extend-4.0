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

namespace Foam
{

void binaryTree::scan()
{
    
    Info << "binaryTree::scan()" << endl;

    if (root_->elementLeft_ !=NULL)
    {
        Info << "root_.elementLeft_" << root_->elementLeft_->v0() << endl;
    }

    if (root_->elementRight_ !=NULL)
    {
        Info << "root_.elementRight_ =" << root_->elementRight_->v0() << endl;
    }
    
    
    forAll(nodeList_,nodei)
    {
        if (nodeList_[nodei]->elementLeft_ !=NULL)
        {
            Info<< "nodeList_[" << nodei
                << "].elementLeft_() ="
                << nodeList_[nodei]->elementLeft_->v0()
                << endl;

            if(nodeList_[nodei]->elementLeft_->node_ ==NULL)
            {
                Info<< "nodeList_[" << nodei
                    << "].elementLeft_() = "
                    << nodeList_[nodei]->elementLeft_->v0()
                    << " doesn't point to any node" << endl;
            }
            
        }
        if (nodeList_[nodei]->elementRight_ !=NULL)
        {
            Info<< "nodeList_[" << nodei
                << "].elementRight_() ="
                << nodeList_[nodei]->elementRight_->v0()
                << endl;

            if(nodeList_[nodei]->elementRight_->node_ ==NULL)
            {
                Info<< "nodeList_[" << nodei
                    << "].elementRight_() ="
                    << nodeList_[nodei]->elementRight_->v0()
                    << " doesn't point to any node" << endl;
            }
        }
        
    }    

    forAll(chemPointList_,pointi)
    {
        if
        (
            chemPointList_[pointi]->node_->elementRight_ == NULL
         && chemPointList_[pointi]->node_->elementLeft_ == NULL
        )
        {
            Info<< "chemPointList_[" << pointi
                << "].node_ "
                << " doesn't point to any node, value = "
                << chemPointList_[pointi]->v0() << endl;
        }
    }    
    Info << "binaryTree::scan() END" << endl;
}


}
