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
    scalar rhostd
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
}


// Construct from components
// CL: needed for operator+
inline mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong
(
    const soaveRedlichKwong& srk,
    label numOfComp,
    DynamicList<scalar> weigths,
    DynamicList<mixtureSoaveRedlichKwong*> mixtureComponents
)
: 
    soaveRedlichKwong(srk),
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
        srk1.numOfComp+1, weigths, mixtureComponents);
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
       srk.mixtureComponents,srk.a0_,srk.b_,srk.Tcrit_,srk.n_,srk.rhostd_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
