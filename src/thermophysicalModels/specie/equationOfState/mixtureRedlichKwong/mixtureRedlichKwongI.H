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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "mixtureRedlichKwong.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
// CL: needed for operator* 
inline mixtureRedlichKwong::mixtureRedlichKwong
(
    const redlichKwong& rK,
    label numOfComp,
    std::vector<scalar> weigths,
    std::vector<mixtureRedlichKwong*> mixtureComponents,
    scalar a,
    scalar b,
    scalar rhostd
)
: 
    redlichKwong(rK),
    numOfComp(numOfComp),
    weigths(weigths),
    mixtureComponents(mixtureComponents)
{
    a_=a;
    b_=b;
    rhostd_=rhostd;
}

// Construct from components
// CL: needed for operator+
inline mixtureRedlichKwong::mixtureRedlichKwong
(
    const redlichKwong& rK,
    label numOfComp,
    std::vector<scalar> weigths,
    std::vector<mixtureRedlichKwong*> mixtureComponents
)
: 
    redlichKwong(rK),
    numOfComp(numOfComp),
    weigths(weigths),
    mixtureComponents(mixtureComponents)    
{
    //CL: update model coefficients
    updateModelCoefficients();
    rhostd_=this->rho(Pstd, Tstd, Pstd/(Tstd*this->R()));
}


// Construct as named copy
inline mixtureRedlichKwong::mixtureRedlichKwong(const word& name, const mixtureRedlichKwong& rK)
:
    redlichKwong(name, rK)
{}


// Construct and return a clone
inline autoPtr<mixtureRedlichKwong> mixtureRedlichKwong::clone() const
{
    return autoPtr<mixtureRedlichKwong>(new mixtureRedlichKwong(*this));
}


// Selector from Istream
inline autoPtr<mixtureRedlichKwong> mixtureRedlichKwong::New(Istream& is)
{
    return autoPtr<mixtureRedlichKwong>(new mixtureRedlichKwong(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//CL: updates the coefficents of the model after the final construction of the mixture
//CL: uses the van der waals mixing rule
inline void mixtureRedlichKwong::updateModelCoefficients() const
{  
    label i,j;
      
    a_=0;
    b_=0;

    for (i=0;i<numOfComp;i++)
    {
        //CL: getting a for the mixture 
        //CL: using van der waals mixing rule
        //CL: TODO: Include possibility to use coupling cofficients k_ij (see paper reference)
        for (j=0;j<numOfComp;j++)
        {
            a_=a_+weigths[i]*weigths[j]*pow(mixtureComponents[i]->a()*mixtureComponents[j]->a(),0.5)*this->W()*this->W();
        }
  
        //CL: getting b for the mixture 
        //CL: using van der waals mixing rule
        b_=b_+weigths[i]*mixtureComponents[i]->b()*this->W();
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


inline mixtureRedlichKwong operator+
(
    const mixtureRedlichKwong& rK1,
    const mixtureRedlichKwong& rK2
)
{
    //CL:save both a pointer to the objector rK1 and rK2 and the weights in the new object
    std::vector<scalar> weigths=rK1.weigths;
    std::vector<mixtureRedlichKwong*> mixtureComponents=rK1.mixtureComponents;

    //CL: Getting the new weigths and mixtureComponents lists, 
    //CL: Saving the object pointer and weigths of rK2 (which is @ mixtureComponents[0] and weigths[0]) at the numOfComp's value of the new object
    weigths.push_back(rK2.weigths[0]);
    mixtureComponents.push_back(rK2.mixtureComponents[0]);

    return mixtureRedlichKwong(static_cast<const redlichKwong&>(rK1)+static_cast<const redlichKwong&>(rK2),
        rK1.numOfComp+1, weigths, mixtureComponents);
}


inline mixtureRedlichKwong operator*
(
    const scalar s,
    const mixtureRedlichKwong& rK
)
{
    //CL: saving the "concentraction" of the component of the mixture in the vector weights
    std::vector<scalar> weigths=rK.weigths;
    weigths.push_back(s*rK.nMoles()); 

   return mixtureRedlichKwong(s*static_cast<const redlichKwong&>(rK), rK.numOfComp, weigths, rK.mixtureComponents,rK.a_,rK.b_,rK.rhostd_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
