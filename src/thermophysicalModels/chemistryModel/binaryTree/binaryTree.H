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

    Data storage of the chemistryOnLineLibrary according to a binary 
    tree structure.
    
          0 (root node)
        /   \
      0      0
    /   \  /   \
   L    R L     0
              /   \
             L     R
             
L: elementLeft_
R: elementRight_
              
\*---------------------------------------------------------------------------*/

#ifndef binary_Tree_H
#define binary_Tree_H

#include "binaryNode.H"      
#include "chemPoint.H"
#include "scalarField.H"
#include "List.H"

namespace Foam
{

class chemistryOnlineLibrary;

class binaryTree
{

private:

        //- Reference to the root node
        binaryNode start_;

        //- Root node of the binary tree
        binaryNode *root_;

        //- Maximum number of elements in the binary tree
        label maxElements_;

        //- Nodes of the tree are also stored in a list
        DynamicList<binaryNode*> nodeList_;

        //- Chempoints of the tree are also stored in a list, used to find the less used element, to replace it when
        //  the library is full 
        DynamicList<chemPoint*> chemPointList_;

        //- Tolerance for the ellipsoid of accuracy  
        scalar tolerance_;
        
        //- Reference to the chemistry on line library
        chemistryOnlineLibrary& onlineLibrary_;

        //- Functions to travel the binary tree: on the left or on the right of a node
        inline bool onLeft(const chemPoint& x, binaryNode *t);
        inline bool onLeft(const scalarField& x, binaryNode *t);
        inline bool onRight(const chemPoint& x, binaryNode *t);
        inline bool onRight(const scalarField& x, binaryNode *t);

        //- Insert a new element on the binary tree
        void insert(const chemPoint& x, binaryNode *t);
        
        //- Insert a new element on the left of the node
        void insertElementLeft(const chemPoint& x, binaryNode *t);

        //- Insert a new element on the right of the node
        void insertElementRight(const chemPoint& x, binaryNode *t);
        
        //- Insert a new node in the tree
        void insertNewNode(const chemPoint& x, binaryNode *t);

        //- Insert a new node on the left of a node
        void insertNewNodeRight(const chemPoint& x, binaryNode *t);

        //- Insert a new node in the right of a node
        void insertNewNodeLeft(const chemPoint& x, binaryNode *t);
        
        //- Check the tree
        void check();
        
        //- Find the less used point in the tree and return a reference to it
        chemPoint& findLessUsed(DynamicList<chemPoint*> pointList);
        chemPoint& findLessUsed(DynamicList<binaryNode*> nodeList);
        chemPoint& findLessUsed(const chemPoint& p, DynamicList<binaryNode*> nodeList);
        
        //- Clean the binary tree
        void clean(chemPoint& oldPoint, binaryNode *oldNode, const chemPoint& newPoint, binaryNode *parentNode);

        //- add a new element to the binary tree, used for replacing
        void add(const chemPoint& newPoint, binaryNode *parentNode, chemPoint& oldPoint, binaryNode *oldNode);

        //- add a new node to the binary tree, used for replacing
        void addNewNode(const chemPoint & x, binaryNode *t, chemPoint& oldPoint, binaryNode *oldNode);

        //- add a new node on the left, used for replacing
        void addNewNodeLeft(const chemPoint & x, binaryNode *t, chemPoint& oldPoint, binaryNode *oldNode);

        //- add a new node on the right, used for replacing
        void addNewNodeRight(const chemPoint & x, binaryNode *t, chemPoint& oldPoint, binaryNode *oldNode);

        //- is the binary tree full?
        bool isListFull();

        //- find the closest element to x
        chemPoint* findClosest(const scalarField& x, binaryNode *t);
        
        //- scale the composition according to the scale factors
        scalarField scaleComposition(const scalarField& x) const;
   
public:

        //- Constructors
        
        //- Construct from dictionary and chemistryOnLineLibrary
        binaryTree
        (
            const dictionary treeDict,
            chemistryOnlineLibrary& onlineLibrary
        );
        
        //- Access
        chemistryOnlineLibrary& onlineLibrary()
        {
            return onlineLibrary_;
        }
        
        //- Edit

        inline void insert(const chemPoint& x);

        //-Search

        inline chemPoint* findClosest(const scalarField& x);

        //- Check
        
        inline label size();        
        void scan();
        
        //- CleanAll
        void cleanAll();

        //- ListFull
        bool ListFull();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "binaryTreeI.H"

#endif
