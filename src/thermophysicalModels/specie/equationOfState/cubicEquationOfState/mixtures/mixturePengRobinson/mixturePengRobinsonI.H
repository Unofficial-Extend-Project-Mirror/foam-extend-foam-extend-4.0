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
    scalar rhostd,
    scalar rhoMin,
    scalar rhoMax
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
    rhoMin_=rhoMin;
    rhoMax_=rhoMax;
}

// Construct from components
// CL: needed for operator+
inline mixturePengRobinson::mixturePengRobinson
(
    const pengRobinson& pr,
    label numOfComp,
    DynamicList<scalar> weigths,
    DynamicList<mixturePengRobinson*> mixtureComponents,
    scalar rhoMin,
    scalar rhoMax
)
: 
    pengRobinson(pr),
    numOfComp(numOfComp),
    weigths(weigths),
    mixtureComponents(mixtureComponents),
    singleComponent(0)
    
{
    TSave=0.0;
    rhoMin_=rhoMin;
    rhoMax_=rhoMax;
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
            label i,j,k;
 
            aSave=0;
            daSave=0;
            d2aSave=0;
            b_=0;

            for (i=0;i<numOfComp;i++)
            { 
                //CL: getting a(T), dadT(T) and d2adT2(T) for the mixture 
                //CL: using van der waals mixing rule
                for (j=0;j<numOfComp;j++)
                {

                    //CL: use mixture correction cofficients k_ij (see paper reference)
                    if(mixtureComponents[0]->realMixtureCorr_==true)
                    {
                        // first and second order temperature derivative of the van der waals mixing rule for a(T) 
                        if(i==j)
                        {
                            aSave=aSave+amix(T,i,j);
                            daSave=daSave+dadTmix(T,i,j);
                            d2aSave=d2aSave+d2adT2mix(T,i,j);
                        }
                        else 
                        {
                            //CL: gives the position of the correction factor in the vector realMixtureCorrCoef_
                            k=i+j-1;

                            aSave=aSave+amix(T,i,j)*(1-mixtureComponents[0]->realMixtureCorrCoef_[k]);
                            daSave=daSave+dadTmix(T,i,j)*(1-mixtureComponents[0]->realMixtureCorrCoef_[k]);
                            d2aSave=d2aSave+d2adT2mix(T,i,j)*(1-mixtureComponents[0]->realMixtureCorrCoef_[k]);  
                        }    

                    }
                    else
                    {
                        aSave=aSave+amix(T,i,j);
                        daSave=daSave+dadTmix(T,i,j);
                        d2aSave=d2aSave+d2adT2mix(T,i,j);
                    }           
                }

                //CL: getting b for the mixture 
                //CL: using van der waals mixing rule
                b_=b_+weigths[i]*mixtureComponents[i]->b()*this->W();

                //CL: saving the temperature at which the mixture coefficients are valid
                TSave=T;
            }
 
            b2_=b_*b_;
            b3_=pow(b_,3);
            b4_=pow(b_,4);
            b5_=pow(b_,5);
        }
    }
}

//CL: returns part of  "a" for the mixture
//CL: Important: a is a function of T
//CL: TODO: RECHECK
inline scalar mixturePengRobinson::amix(const scalar T, const label i, const label j) const
{

    return weigths[i]*weigths[j]*this->W()*this->W()*
               pow(mixtureComponents[i]->a(T)*mixtureComponents[j]->a(T),0.5);
}

//CL: returns part of the first derivative of "a" for the mixture
//CL: Important: a is a function of T
//CL: TODO: RECHECK
inline scalar mixturePengRobinson::dadTmix(const scalar T, const label i, const label j) const
{

    return 0.5*weigths[i]*weigths[j]*this->W()*this->W()
               *(pow(mixtureComponents[i]->a(T)/mixtureComponents[j]->a(T),0.5)*mixtureComponents[j]->dadT(T)
                   +pow(mixtureComponents[j]->a(T)/mixtureComponents[i]->a(T),0.5)*mixtureComponents[i]->dadT(T));
}

//CL: returns part of the second derivative of "a" for the mixture
//CL: Important: a is a function of T
//CL: TODO: RECHECK
inline scalar mixturePengRobinson::d2adT2mix(const scalar T, const label i, const label j) const
{

    return  0.5*weigths[i]*weigths[j]*this->W()*this->W()*
            (
                pow(mixtureComponents[i]->a(T)/mixtureComponents[j]->a(T),0.5)*mixtureComponents[j]->d2adT2(T)
                    +pow(mixtureComponents[j]->a(T)/mixtureComponents[i]->a(T),0.5)*mixtureComponents[i]->d2adT2(T)
                        +pow(mixtureComponents[i]->a(T)*mixtureComponents[j]->a(T),-0.5)
                            *mixtureComponents[i]->dadT(T)*mixtureComponents[j]->dadT(T)
                                -0.5*pow(mixtureComponents[i]->a(T)/mixtureComponents[j]->a(T),0.5)*1/mixtureComponents[j]->a(T)
                                    *pow(mixtureComponents[j]->dadT(T),2)
                                        -0.5*pow(mixtureComponents[j]->a(T)/mixtureComponents[i]->a(T),0.5)
                                            *1/mixtureComponents[i]->a(T)*pow(mixtureComponents[i]->dadT(T),2)
            );
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

    //CL: setting TSave=0--> results that all coefficients (a, da, d2a) are calculated in updateModelCoefficients
    TSave=0.0;
    singleComponent=0;

    //CL:setting rho boundaries
    rhoMin_=min(rhoMin_,pr.rhoMin_);
    rhoMax_=max(rhoMax_,pr.rhoMax_);

    //CL: calculating rho @ std
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
        pr1.numOfComp+1, weigths, mixtureComponents, min(pr1.rhoMin_,pr2.rhoMin_), max(pr1.rhoMax_,pr2.rhoMax_));
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
       pr.mixtureComponents,pr.a0_,pr.b_,pr.Tcrit_,pr.n_,pr.rhostd_, pr.rhoMin_, pr.rhoMax_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
