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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "chemPoint.H"
#include "binaryNode.H"

// #define debugEOA
// #define debugGrow
// #define debugCheckSolution

namespace Foam
{

chemPoint::chemPoint
(
    const scalarField& v0,
    const scalarField& r,
    const scalarField& tolerances,
    const scalarField& tolerancesSolutions,
    const scalar& absErr,
    const Switch& logT,
    const scalar& deltaT
)
:
    node_(NULL),
    nUsed_(0),
    v0_(v0),
    r_(r),
    EOA_(tolerances.size()),
    solutionsEOA_(tolerancesSolutions.size()),
    absErr_(absErr),
    logT_(logT),
    deltaT_(deltaT),
    timeEOA_(0)
    
{


    forAll(tolerances, i)
    {
        if(tolerances[i] == 0)
        {
            EOA_[i] = GREAT;
        }
        else
        {
            EOA_[i] = absErr * tolerances[i];        
        }
    }
    
    forAll(tolerancesSolutions, i)
    {
        if(tolerances[i] == 0)
        {
            solutionsEOA_[i] = GREAT;
        }
        else
        {
            solutionsEOA_[i] = absErr * tolerancesSolutions[i];        
        }
    }
    
}


chemPoint::chemPoint
(
    const scalarField& v0,
    const scalarField& r,
    const scalarField& tolerances,
    const scalarField& tolerancesSolutions,
    const scalar& absErr,
    const Switch& logT,
    const scalar& deltaT,
    binaryNode* node
)
:
    node_(node),
    nUsed_(0),
    v0_(v0),
    r_(r),
    EOA_(tolerances.size()),
    solutionsEOA_(tolerancesSolutions.size()),
    absErr_(absErr),
    logT_(logT),
    deltaT_(deltaT),
    timeEOA_(0)

{
    forAll(tolerances, i)
    {
        if(tolerances[i] == 0)
        {
            EOA_[i] = 0.0;
        }
        else
        {
            EOA_[i] = absErr * tolerances[i];        
        }
    }    
    
    forAll(tolerancesSolutions, i)
    {
        if(tolerancesSolutions[i] == 0)
        {
            solutionsEOA_[i] = 0.0;
        }
        else
        {
            solutionsEOA_[i] = absErr * tolerancesSolutions[i];        
        }
    }


}

chemPoint::chemPoint
(
    const chemPoint& p,
    binaryNode* node
)
:
    node_(node),
    nUsed_(0),
    v0_(p.v0()),
    r_(p.r()),
    EOA_(p.EOA()),
    solutionsEOA_(p.solutionsEOA()),
    absErr_(p.absErr()),
    logT_(p.logT()),
    deltaT_(p.deltaT()),
    timeEOA_(p.timeEOA())
{}    

chemPoint::chemPoint
(
    const chemPoint& p
)
:
    node_(NULL),
    nUsed_(0),
    v0_(p.v0()),
    r_(p.r()),
    EOA_(p.EOA()),
    solutionsEOA_(p.solutionsEOA()),
    absErr_(p.absErr()),
    logT_(p.logT()),
    deltaT_(p.deltaT()),
    timeEOA_(p.timeEOA())
{}    

bool chemPoint::inEOA
(
    const scalarField& point,
    const scalarField& Wi,
    const scalar& rhoi, 
    const scalar& deltaT,
    const scalarField& scaleFactor 
)
{
    
    if(nUsed_ < INT_MAX)
    {
        nUsed_++;
    }

    label size = point.size();    
        
    scalar eps2 = 0.0;

    for(label i=0; i < size; i++)
    {

//      Temperature    
        if(i == size-2 && logT_)
        {
            scalar diffAxis = sqr(Foam::scalar(mag(log(point[i])-log(v0_[i]))))/sqr(EOA_[i]) ;
            eps2 += diffAxis;
        }
        else if(i == size-2 && !logT_)
        {
            scalar diffAxis = sqr(Foam::scalar(mag((point[i])-(v0_[i])))) / sqr(EOA_[i]) ;
            eps2 += diffAxis;                        
        }
        if(i == size-1)
        {
            scalar diffAxis = sqr(Foam::scalar(mag((point[i]-v0_[i])))) / sqr(EOA_[i]) ;
            eps2 += diffAxis;
        }
//      Species                
        else
        {
            scalar diffAxis = sqr(Foam::scalar(mag((point[i]-v0_[i])))) / sqr(EOA_[i]) ;
            eps2 += diffAxis;
        }
    }
       
    
    if(eps2 > 1.0)
    {
        return false;
    }
    else
    {            
        return true;
    }
    
}

void chemPoint::grow
(
    const scalarField& v0,
    const scalarField& Wi,
    const scalar& rhoi,
    const scalarField& v, 
    const scalar& deltaT
)
{

    label size = v.size();
    
    for(label i=0; i < size-2; i++)
    {
        if(i == size-2 && logT_)
        {
            scalar diffAxis = Foam::scalar(mag(log(v0[i])-log(v0_[i]))) - EOA_[i] ;

            if(diffAxis > Foam::VSMALL)
            {
                EOA_[i] = mag(log(v0[i])-log(v0_[i]));  
            }
        }
        else if(i == size-2 && !logT_)
        {
            scalar diffAxis = Foam::scalar(mag((v0[i])-(v0_[i]))) - EOA_[i] ;

            if(diffAxis > Foam::VSMALL)
            {
                EOA_[i] = mag((v0[i])-(v0_[i]));  
            }
        }
        else
        {
            scalar diffAxis = Foam::scalar(mag((v0[i]-v0_[i]))) - EOA_[i] ;

            if(diffAxis > Foam::VSMALL)
            {
                EOA_[i] = mag((v0[i]-v0_[i]));  
            }
            
        }
    }
    
    timeEOA_ = Foam::scalar(mag(deltaT_ - deltaT));
        
}

bool chemPoint::checkSolution
(
    const scalarField& v0,
    const scalarField& v,
    const scalarField& Wi,
    const scalar& rhoi, 
    const scalar& T,
    const scalar& p,
    const scalar& deltaT,
    const scalarField& tolerances
)
{

    label size = v.size();
    
    scalarField s = v;

    s.setSize(size+2);

//    s[size-2] = T;
//    s[size-1] = p;

    s[s.size()-2] = T;
    s[s.size()-1] = p;
    size = s.size();
    
    scalar eps2 = 0.0;


    for(label i=0; i < size; i++)
    {
        if(i == size-2 && logT_)
        {
            scalar diffAxis = (Foam::scalar(mag(log(s[i])-log(r_[i])))) / (solutionsEOA_[i]) ;
            eps2 += sqr(diffAxis);            
        }
        else if(i == size-2 && !logT_)
        {
            scalar diffAxis = (Foam::scalar(mag(s[i]-r_[i]))) / (solutionsEOA_[i]) ;
            eps2 += sqr(diffAxis);                        
        }        
        else if(i == size-1)
        {
            scalar diffAxis = (Foam::scalar(mag(s[i]-r_[i]))) / (solutionsEOA_[i]) ;
            eps2 += sqr(diffAxis);                        
        }        
        else
        {
            scalar diffAxis = (Foam::scalar(mag((s[i]*Wi[i]/rhoi-r_[i])))) / (solutionsEOA_[i]) ;
            eps2 += sqr(diffAxis);
        }
    }

    if(eps2 > 1.0)
    {
        return false;
    }
    else
    {
        // if the solution is in the ellipsoid of accuracy, grow it!
        grow(v0, Wi, rhoi, s, deltaT);
        return true;    
    }


}

void chemPoint::setFree()
{
    node_ = NULL;
    nUsed_ = 0;
    /*
        v0_.clear();        
        r_.clear();
        EOA_.clear();        
        solutionsEOA_.clear();
        absErr_ = 0;
        logT_ = 0;
        deltaT_ = 0;
        timeEOA_ = 0; 
    */        
}

void chemPoint::clearData()
{
        nUsed_ = 0;
//        Info << "nUsed" << endl;       
        v0_.clear();        
//        Info << "v0.clear" << endl;       
        r_.clear();
//        Info << "r_.clear()" << endl;
        EOA_.clear();        
//        Info << "EOA_.clear()" << endl;
        solutionsEOA_.clear();
//        Info << "solutionsEOA_.clear()" << endl;
        absErr_ = 0;
//        Info << "absErr_ = 0" << endl;
//        logT_ = 0;
//        Info << "logT_ = 0" << endl;
        deltaT_ = 0;
//        Info << "deltaT_ = 0;" << endl;
        timeEOA_ = 0; 
//        node_ = NULL;
//        Info << "chemPoint::clearData():: ENDDDDDD!" << endl;

}


}
