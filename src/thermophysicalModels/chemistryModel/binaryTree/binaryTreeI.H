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

namespace Foam
{

inline void binaryTree::insert
( 
    const chemPoint& x 
) 
{
    insert(x, root_);
}

inline chemPoint* binaryTree::findClosest
( 
    const scalarField& x 
) 
{
    return findClosest(x, root_);
}


inline bool binaryTree::onLeft
(
    const chemPoint& x, 
    binaryNode *t
)
{
    const scalarField& v = scaleComposition(x.v0());
        
    if(sum(v*t->vT()) < t->a())
    {
        // element on the left side
        return true;  
    }
    else
    {
        return false;
    }
}


inline bool binaryTree::onRight
(
    const chemPoint& x, 
    binaryNode *t
)
{
    const scalarField& v = scaleComposition(x.v0());

    if(sum(v*t->vT()) > t->a())
    {
        // element on the right side
        return true;  
    }
    else
    {
        return false;
    }
}


inline bool binaryTree::onLeft
(
    const scalarField& x, 
    binaryNode *t
)
{
    scalarField xS = scaleComposition(x);

    if(sum(xS*t->vT()) < t->a())
    {
        // element on the left side
        return true;  
    }
    else
    {
        return false;
    }
}

inline bool binaryTree::onRight
(
    const scalarField& x, 
    binaryNode *t
)
{
    scalarField xS = scaleComposition(x);

    if(sum(xS*t->vT()) > t->a())
    {
        // element on the right side
        return true;  

    }
    else
    {
        return false;
    }
}

inline label binaryTree::size()
{
    return chemPointList_.size();
}


}
