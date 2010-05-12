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
    
    Node of the binary tree

\*---------------------------------------------------------------------------*/

#ifndef BINARY_NODE_H
#define BINARY_NODE_H

#include "chemPoint.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class binaryNode Declaration
\*---------------------------------------------------------------------------*/

class binaryNode
{

public:

        //- Element on the left
        chemPoint* elementLeft_;

        //- Element on the right
        chemPoint* elementRight_;

        //- Node which follows on the left
        binaryNode* left_;

        //- Node which follows on the right
        binaryNode* right_;

        //- Upper node
        binaryNode* down_;


private:
        
        //- vT and a are used for travelling the node during
        //  the search algorithm
        scalarField vT_;
        scalar a_;
        
        //- vLeft_ and vRight_ are stored in the node, in order to
        //  update vT and a when a new element is inserted
        scalarField vLeft_;
        scalarField vRight_;
            
public:

    // Constructors
        
        //- Construct null
        binaryNode();
        
        //- Construct from components
        binaryNode
        (
            chemPoint* elementLeft,
            chemPoint* elementRight,
            binaryNode* left,
            binaryNode* right,
            binaryNode* down
        );
        
        //- Construct from another binary node
        binaryNode
        (
            binaryNode *bn
        );
        

    // Member functions

        //- Access

        inline chemPoint* elementLeft()
        {
            return elementLeft_;
        }

        inline chemPoint* elementRight()
        {
            return elementRight_;
        }

        inline binaryNode* left()
        {
            return left_;
        }

        inline binaryNode* right()
        {
            return right_;
        }

        inline binaryNode* down()
        {
            return down_;
        }

        inline scalarField& vLeft()
        {
            return vLeft_;
        }

        inline scalarField& vRight()
        {
            return vRight_;
        }


        //- Edit

        //- Set the binary node free from element left and element right
        void setFree();

        //- Topology

        inline scalarField& vT()
        {
            return vT_;
        }

        inline scalar& a()
        {
            return a_;        
        }

        //- Update vT and a in the binary node
        void updateNode();

        void clearData();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
