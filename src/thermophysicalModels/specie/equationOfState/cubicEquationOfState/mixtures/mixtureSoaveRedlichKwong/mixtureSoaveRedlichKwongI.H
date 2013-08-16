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

#include "mixtureSoaveRedlichKwong.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
// CL: needed for operator* 
inline mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong
(
    const soaveRedlichKwong& srk,
    label numOfComp,
    DynamicList<scalar> weigths,
    DynamicList<mixtureSoaveRedlichKwong*> mixtureComponents,
    scalar a0,
    scalar b,
    scalar Tcrit,
    scalar n,
    scalar rhostd,
    scalar rhoMin,
    scalar rhoMax
)
: 
    soaveRedlichKwong(srk),
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
inline mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong
(
    const soaveRedlichKwong& srk,
    label numOfComp,
    DynamicList<scalar> weigths,
    DynamicList<mixtureSoaveRedlichKwong*> mixtureComponents,
    scalar rhoMin,
    scalar rhoMax
)
: 
    soaveRedlichKwong(srk),
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
inline mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong(const word& name, const mixtureSoaveRedlichKwong& srk)
:
    soaveRedlichKwong(name, srk)
{}


// Construct and return a clone
inline autoPtr<mixtureSoaveRedlichKwong> mixtureSoaveRedlichKwong::clone() const
{
    return autoPtr<mixtureSoaveRedlichKwong>(new mixtureSoaveRedlichKwong(*this));
}


// Selector from Istream
inline autoPtr<mixtureSoaveRedlichKwong> mixtureSoaveRedlichKwong::New(Istream& is)
{
    return autoPtr<mixtureSoaveRedlichKwong>(new mixtureSoaveRedlichKwong(is));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//CL: updates the coefficents of the model after the final construction of the mixture
//CL: uses the van der waals mixing rule
inline void mixtureSoaveRedlichKwong::updateModelCoefficients(const  scalar T) const
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
            b5_=pow(b_,5);
        }
    }
}

//CL: returns part of  "a" for the mixture
//CL: Important: a is a function of T
//CL: TODO: RECHECK
inline scalar mixtureSoaveRedlichKwong::amix(const scalar T, const label i, const label j) const
{

    return weigths[i]*weigths[j]*this->W()*this->W()*
               pow(mixtureComponents[i]->a(T)*mixtureComponents[j]->a(T),0.5);
}

//CL: returns part of the first derivative of "a" for the mixture
//CL: Important: a is a function of T
//CL: TODO: RECHECK
inline scalar mixtureSoaveRedlichKwong::dadTmix(const scalar T, const label i, const label j) const
{

    return 0.5*weigths[i]*weigths[j]*this->W()*this->W()
               *(pow(mixtureComponents[i]->a(T)/mixtureComponents[j]->a(T),0.5)*mixtureComponents[j]->dadT(T)
                   +pow(mixtureComponents[j]->a(T)/mixtureComponents[i]->a(T),0.5)*mixtureComponents[i]->dadT(T));
}

//CL: returns part of the second derivative of "a" for the mixture
//CL: Important: a is a function of T
//CL: TODO: RECHECK
inline scalar mixtureSoaveRedlichKwong::d2adT2mix(const scalar T, const label i, const label j) const
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
inline scalar mixtureSoaveRedlichKwong::p(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::p(rho,T);
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar mixtureSoaveRedlichKwong::dpdv(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::dpdv(rho,T);
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar mixtureSoaveRedlichKwong::dpdT(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::dpdT(rho,T);
}


//Real  deviative dv/dT at constant pressure 
//(molar values)
inline scalar mixtureSoaveRedlichKwong::dvdT(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::dvdT(rho,T);
}


//Real  deviative dv/dp at constant temperature 
//(molar values)
inline scalar mixtureSoaveRedlichKwong::dvdp(const scalar rho, const scalar T)  const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::dvdp(rho,T);
}


//needed to calculate the internal energy
//(molar values)
inline scalar mixtureSoaveRedlichKwong::integral_p_dv
(
    const scalar rho,
    const scalar T
) const
{ 
    updateModelCoefficients(T);
    return soaveRedlichKwong::integral_p_dv(rho,T);
}


//needed to calculate the entropy
//(molar values)
inline scalar mixtureSoaveRedlichKwong::integral_dpdT_dv
(
    const scalar rho,
    const scalar T
) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::integral_dpdT_dv(rho,T);
}


//(molar values)
inline scalar mixtureSoaveRedlichKwong::d2pdT2(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::d2pdT2(rho,T);
} 


//(molar values)
inline scalar mixtureSoaveRedlichKwong::d2pdv2(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::d2pdv2(rho,T);
}


//(molar values)
inline scalar mixtureSoaveRedlichKwong::d2vdT2(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::d2vdT2(rho,T);
}


//(molar values)
inline scalar mixtureSoaveRedlichKwong::d2pdvdT(const scalar rho, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::d2pdvdT(rho,T);
}


// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar mixtureSoaveRedlichKwong::integral_d2pdT2_dv
(
    const scalar rho,
    const scalar T
)  const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::integral_d2pdT2_dv(rho,T);
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar mixtureSoaveRedlichKwong::isobarExpCoef
(
    const scalar rho,
    const scalar T
) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::isobarExpCoef(rho,T);
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar mixtureSoaveRedlichKwong::isothermalCompressiblity
(
    const scalar rho,
    const scalar T
) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::isothermalCompressiblity(rho,T);
}


//- Return density [kg/m^3]on 
inline scalar mixtureSoaveRedlichKwong::rho
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    updateModelCoefficients(T);

    return soaveRedlichKwong::rho(p,T,rho0);
}


//- Return density [kg/m^3]on 
inline scalar mixtureSoaveRedlichKwong::rho(const scalar p, const scalar T) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::rho(p,T);
}


//- Return compressibility drho/dp at T=constant [s^2/m^2]
inline scalar mixtureSoaveRedlichKwong::psi(const scalar rho, const scalar T) const
{
    this->updateModelCoefficients(T);
    return soaveRedlichKwong::psi(rho,T);
}


//- Return compression factor []
inline scalar mixtureSoaveRedlichKwong::Z
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    updateModelCoefficients(T);
    return soaveRedlichKwong::Z(p,T,rho0);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


inline void mixtureSoaveRedlichKwong::operator+=(const mixtureSoaveRedlichKwong& srk)
{
    //CL: Resizes the DynamicList if number of compontens is larger than the DynamicList
    if((weigths.size()<numOfComp+1)||(mixtureComponents.size()<numOfComp+1))
    {
        weigths.setSize(2*numOfComp);
        mixtureComponents.setSize(2*numOfComp);
    }

    //Cl: oldClass+=srk 
    //CL: Saving the object pointer and weigths of pr (which is @ mixtureComponents[0] and weigths[0]) at the numOfComp's value of the oldClass 
    weigths[numOfComp]=srk.weigths[0];
    mixtureComponents[numOfComp]=srk.mixtureComponents[0];

    soaveRedlichKwong::operator+=(srk);

    //CL: increase number of Components by 1
    numOfComp=numOfComp+1;

    //CL: setting TSave=0--> results that all coefficients (a, da, d2a) are calculated in updateModelCoefficients
    TSave=0.0;
    singleComponent=0;

    //CL:setting rho boundaries
    rhoMin_=min(rhoMin_,srk.rhoMin_);
    rhoMax_=max(rhoMax_,srk.rhoMax_);

    //CL: calculating rho @ std
    rhostd_=this->rho(Pstd, Tstd, Pstd*this->W()/(Tstd*this->R()));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


inline mixtureSoaveRedlichKwong operator+
(
    const mixtureSoaveRedlichKwong& srk1,
    const mixtureSoaveRedlichKwong& srk2
)
{
    //CL:save both a pointer to the objector srk1 and srk2 and the weights in the new object
    DynamicList<scalar> weigths=srk1.weigths;
    DynamicList<mixtureSoaveRedlichKwong*> mixtureComponents=srk1.mixtureComponents;

    //CL: Resizes the DynamicList if number of compontens is larger than the DynamicList
    if((weigths.size()<srk1.numOfComp+1)||(mixtureComponents.size()<srk1.numOfComp+1))
    {
        weigths.setSize(2*srk1.numOfComp);
        mixtureComponents.setSize(2*srk1.numOfComp);
    }

    //CL: Getting the new weigths and mixtureComponents lists, 
    //CL: Saving the object pointer and weigths of srk2 (which is @ mixtureComponents[0] and weigths[0]) at the numOfComp's value of the new object
    weigths[srk1.numOfComp]=srk2.weigths[0];
    mixtureComponents[srk1.numOfComp]=srk2.mixtureComponents[0];

    return mixtureSoaveRedlichKwong(static_cast<const soaveRedlichKwong&>(srk1)+static_cast<const soaveRedlichKwong&>(srk2),
        srk1.numOfComp+1, weigths, mixtureComponents, min(srk1.rhoMin_,srk2.rhoMin_), max(srk1.rhoMax_,srk2.rhoMax_));
}


inline mixtureSoaveRedlichKwong operator*
(
    const scalar s,
    const mixtureSoaveRedlichKwong& srk
)
{
    //CL: saving the "concentraction" of the component of the mixture in the vector weights
    //CL: saved at the Position "numOfComp-1"
    DynamicList<scalar> weigths=srk.weigths;
    weigths[srk.numOfComp-1]=s*srk.nMoles(); 

   return mixtureSoaveRedlichKwong(s*static_cast<const soaveRedlichKwong&>(srk), srk.numOfComp, weigths,
       srk.mixtureComponents,srk.a0_,srk.b_,srk.Tcrit_,srk.n_,srk.rhostd_, srk.rhoMin_, srk.rhoMax_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
