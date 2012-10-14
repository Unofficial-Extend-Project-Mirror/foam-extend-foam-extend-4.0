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

#include "mixturePengRobinson.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
// CL: needed for operator* 
inline mixturePengRobinson::mixturePengRobinson
(
    const pengRobinson& pr,
    label numOfComp,
    DynamicList<scalar> weigths,
    DynamicList<mixturePengRobinson*> mixtureComponents,
    scalar a0,
    scalar b,
    scalar Tcrit,
    scalar n,
    scalar rhostd
)
: 
    pengRobinson(pr),
    numOfComp(numOfComp),
    weigths(weigths),
    mixtureComponents(mixtureComponents),
    singleComponent(1)
{
    a0_=a0;
    b_=b;
    Tcrit_=Tcrit; 
    n_=n;
    rhostd_=rhostd;
}

// Construct from components
// CL: needed for operator+
inline mixturePengRobinson::mixturePengRobinson
(
    const pengRobinson& pr,
    label numOfComp,
    DynamicList<scalar> weigths,
    DynamicList<mixturePengRobinson*> mixtureComponents
)
: 
    pengRobinson(pr),
    numOfComp(numOfComp),
    weigths(weigths),
    mixtureComponents(mixtureComponents),
    singleComponent(0)
    
{
    TSave=0.0;
    //CL: following three DynamicList have no size in the new object
    //CL: resize DynamicList
    aComponents.resize(numOfComp);
    daComponents.resize(numOfComp);
    d2aComponents.resize(numOfComp);
    rhostd_=this->rho(Pstd, Tstd, Pstd*this->W()/(Tstd*this->R()));
}


// Construct as named copy
inline mixturePengRobinson::mixturePengRobinson(const word& name, const mixturePengRobinson& pr)
:
    pengRobinson(name, pr)
{}


// Construct and return a clone
inline autoPtr<mixturePengRobinson> mixturePengRobinson::clone() const
{
    return autoPtr<mixturePengRobinson>(new mixturePengRobinson(*this));
}


// Selector from Istream
inline autoPtr<mixturePengRobinson> mixturePengRobinson::New(Istream& is)
{
    return autoPtr<mixturePengRobinson>(new mixturePengRobinson(is));
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//CL: updates the coefficents of the model after the construction of the mixture
//CL: uses the van der waals mixing rule
inline void mixturePengRobinson::updateModelCoefficients(const  scalar T) const
{  
    // CL: to make sure that the coefficents are only updated if the mixture has more than 1 component 
    if (singleComponent==0)
    {    
        // Checking if the mixture coefficient were already calculated for this temperature
        if(TSave!=T)
        {
            label i,j;
 
            aSave=0;
            daSave=0;
            d2aSave=0;
            b_=0;

            //CL: filling vector a, dadT and d2adT2
            for (i=0;i<numOfComp;i++)
            {
                aComponents[i]=mixtureComponents[i]->a(T);
                daComponents[i]=mixtureComponents[i]->dadT(T);
                d2aComponents[i]=mixtureComponents[i]->d2adT2(T);
            }

            for (i=0;i<numOfComp;i++)
            { 
                //CL: getting a(T), dadT(T) and d2adT2(T) for the mixture 
                //CL: using van der waals mixing rule
                //CL: TODO: Include possibility to use coupling cofficients k_ij (see paper reference)
                for (j=0;j<numOfComp;j++)
                {
                    aSave=aSave+weigths[i]*weigths[j]*pow(aComponents[i]*aComponents[j],0.5)*this->W()*this->W();
   
                    // first and second order temperature derivative of the van der waals mixing rule for a(T) 
                    if(i==j)
                    {
                        daSave=daSave+weigths[i]*weigths[i]*daComponents[i]*this->W()*this->W();
                        d2aSave=d2aSave+weigths[i]*weigths[i]*d2aComponents[i]*this->W()*this->W();
                    }
                    else 
                    {
                        daSave=daSave+0.5*weigths[i]*weigths[j]*this->W()*this->W()
                            *(pow(aComponents[i]/aComponents[j],0.5)*daComponents[j]+pow(aComponents[j]/aComponents[i],0.5)*daComponents[i]);
                   
                        d2aSave=d2aSave+0.5*weigths[i]*weigths[j]*this->W()*this->W()*
                             (
                                 pow(aComponents[i]/aComponents[j],0.5)*d2aComponents[i]+pow(aComponents[j]/aComponents[i],0.5)*d2aComponents[j]
                                     +0.5*pow(aComponents[i]*aComponents[j],-0.5)*(pow(daComponents[i],2)+pow(daComponents[j],2))
                                         -0.5*daComponents[i]*daComponents[j]*(pow(aComponents[i]/aComponents[j],0.5)/aComponents[j]
                                             +pow(aComponents[j]/aComponents[i],0.5)/aComponents[i])
                             );   
                    }              
                }

                //CL: getting b for the mixture 
                //CL: using van der waals mixing rule
                b_=b_+weigths[i]*mixtureComponents[i]->b()*this->W();

                //CL: saving the temperature at which the mixture coefficients are valid
                TSave=T;
            }
        }
    }

}


//returns the pressure for a given density and temperature    
inline scalar mixturePengRobinson::p(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::p(rho,T);
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar mixturePengRobinson::dpdv(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::dpdv(rho,T);
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar mixturePengRobinson::dpdT(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::dpdT(rho,T);
}


//Real  deviative dv/dT at constant pressure 
//(molar values)
inline scalar mixturePengRobinson::dvdT(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::dvdT(rho,T);
}


//Real  deviative dv/dp at constant temperature 
//(molar values)
inline scalar mixturePengRobinson::dvdp(const scalar rho, const scalar T)  const
{
    updateModelCoefficients(T);
    return pengRobinson::dvdp(rho,T);
}


//needed to calculate the internal energy
//(molar values)
inline scalar mixturePengRobinson::integral_p_dv
(
    const scalar rho,
    const scalar T
) const
{ 
    updateModelCoefficients(T);
    return pengRobinson::integral_p_dv(rho,T);
}


//needed to calculate the entropy
//(molar values)
inline scalar mixturePengRobinson::integral_dpdT_dv
(
    const scalar rho,
    const scalar T
) const
{
    updateModelCoefficients(T);
    return pengRobinson::integral_dpdT_dv(rho,T);
}


//(molar values)
inline scalar mixturePengRobinson::d2pdT2(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::d2pdT2(rho,T);
} 


//(molar values)
inline scalar mixturePengRobinson::d2pdv2(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::d2pdv2(rho,T);
}


//(molar values)
inline scalar mixturePengRobinson::d2vdT2(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::d2vdT2(rho,T);
}


//(molar values)
inline scalar mixturePengRobinson::d2pdvdT(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::d2pdvdT(rho,T);
}


// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar mixturePengRobinson::integral_d2pdT2_dv
(
    const scalar rho,
    const scalar T
)  const
{
    updateModelCoefficients(T);
    return pengRobinson::integral_d2pdT2_dv(rho,T);
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar mixturePengRobinson::isobarExpCoef
(
    const scalar rho,
    const scalar T
) const
{
    updateModelCoefficients(T);
    return pengRobinson::isobarExpCoef(rho,T);
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar mixturePengRobinson::isothermalCompressiblity
(
    const scalar rho,
    const scalar T
) const
{
    updateModelCoefficients(T);
    return pengRobinson::isothermalCompressiblity(rho,T);
}


//- Return density [kg/m^3]on 
inline scalar mixturePengRobinson::rho
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    updateModelCoefficients(T);
    return pengRobinson::rho(p,T,rho0);
}


//- Return density [kg/m^3]on 
inline scalar mixturePengRobinson::rho(const scalar p, const scalar T) const
{
    updateModelCoefficients(T);
    return pengRobinson::rho(p,T);
}


//- Return compressibility drho/dp at T=constant [s^2/m^2]
inline scalar mixturePengRobinson::psi(const scalar rho, const scalar T) const
{
    this->updateModelCoefficients(T);
    return pengRobinson::psi(rho,T);
}


//- Return compression factor []
inline scalar mixturePengRobinson::Z
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    updateModelCoefficients(T);
    return pengRobinson::Z(p,T,rho0);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


inline void mixturePengRobinson::operator+=(const mixturePengRobinson& pr)
{
    //CL: Resizes the DynamicList if number of compontens is larger than the DynamicList
    if((weigths.size()<numOfComp+1)||(mixtureComponents.size()<numOfComp+1))
    {
        weigths.setSize(2*numOfComp);
        mixtureComponents.setSize(2*numOfComp);
    }
  
    //Cl: oldClass+=pr 
    //CL: Saving the object pointer and weigths of pr (which is @ mixtureComponents[0] and weigths[0]) at the numOfComp's value of the oldClass 
    weigths[numOfComp]=pr.weigths[0];
    mixtureComponents[numOfComp]=pr.mixtureComponents[0];

    pengRobinson::operator+=(pr);

    //CL: increase number of Components by 1
    numOfComp=numOfComp+1;

    //CL: resize DynamicList so that they have a defined size in the mixture object
    aComponents.resize(numOfComp);
    daComponents.resize(numOfComp);
    d2aComponents.resize(numOfComp);

    //CL: setting TSave=0--> results that all coefficients (a, da, d2a) are calculated in updateModelCoefficients
    TSave=0.0;
    singleComponent=0;
    rhostd_=this->rho(Pstd, Tstd, Pstd*this->W()/(Tstd*this->R()));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline mixturePengRobinson operator+
(
    const mixturePengRobinson& pr1,
    const mixturePengRobinson& pr2
)
{
    //CL:save both a pointer to the objector pr1 and pr2 and the weights in the new object
    DynamicList<scalar> weigths=pr1.weigths;
    DynamicList<mixturePengRobinson*> mixtureComponents=pr1.mixtureComponents;

    //CL: Resizes the DynamicList if number of compontens is larger than the DynamicList
    if((weigths.size()<pr1.numOfComp+1)||(mixtureComponents.size()<pr1.numOfComp+1))
    {
        weigths.setSize(2*pr1.numOfComp);
        mixtureComponents.setSize(2*pr1.numOfComp);
    }

    //CL: Getting the new weigths and mixtureComponents lists, 
    //CL: Saving the object pointer and weigths of pr2 (which is @ mixtureComponents[0] and weigths[0]) at the numOfComp's value of the new object
    weigths[pr1.numOfComp]=pr2.weigths[0];
    mixtureComponents[pr1.numOfComp]=pr2.mixtureComponents[0];

    return mixturePengRobinson(static_cast<const pengRobinson&>(pr1)+static_cast<const pengRobinson&>(pr2),
        pr1.numOfComp+1, weigths, mixtureComponents);
}


inline mixturePengRobinson operator*
(
    const scalar s,
    const mixturePengRobinson& pr
)
{
    //CL: saving the "concentraction" of the component of the mixture in the vector weights
    //CL: saved at the Position "numOfComp-1"
    DynamicList<scalar> weigths=pr.weigths;
    weigths[pr.numOfComp-1]=s*pr.nMoles(); 

   return mixturePengRobinson(s*static_cast<const pengRobinson&>(pr), pr.numOfComp, weigths,
       pr.mixtureComponents,pr.a0_,pr.b_,pr.Tcrit_,pr.n_,pr.rhostd_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
