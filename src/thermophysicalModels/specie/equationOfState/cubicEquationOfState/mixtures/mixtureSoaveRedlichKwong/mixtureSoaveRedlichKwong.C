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

Description
    Soave Redlich Kwong equation of state for mixtures.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "mixtureSoaveRedlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong(Istream& is)
:
    soaveRedlichKwong(is),
    numOfComp(1),
    singleComponent(1),
    //CL: no real gas mixture correction when stream constructor is used
    realMixtureCorr_(false)
{ 
    //CL: set size of weigths, mixtureComponents ... to 10,
    //CL: when more mixture componentents are used
    //CL: size of the DynamicLis increases automatically
    weigths.setSize(10);
    mixtureComponents.setSize(10);

    //Save a pointer of this object in the mixtureComponents array
    mixtureComponents[0]=this;
    is.check("mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong(Istream& is)");
}

//CL: Constructed needed in OpenFOAM 2.x.x
//CL: Code works fine, but compiling problem in OpenFOAM 1.6.ext
//CL:  because specie has no constructor using dict
/*
mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong(const dictionary& dict)
:
    soaveRedlichKwong(dict),
    numOfComp(1),
    singleComponent(1),
    realMixtureCorr_(dict.subDict("equationOfState").lookupOrDefault("realMixtureCorrection",false))
{ 
    //CL: set size of weigths, mixtureComponents ... to 10,
    //CL: when more mixture componentents are used
    //CL: size of the DynamicLis increases automatically
    weigths.setSize(10);
    mixtureComponents.setSize(10);

    //Save a pointer of this object in the mixtureComponents array
    mixtureComponents[0]=this;

    if(realMixtureCorr_==true)
    {
        //CL: reads number of mixture components
        //CL: this number is needed to calculate the number of correction coefficients
        nCom_=dict.subDict("equationOfState").lookupOrDefault("numberOfComponents",1);

            if (nCom_<1||nCom_==1)
            {
                FatalErrorIn
                (
                    "mixturePengRobinson::mixturePengRobinson(const dictionary& dict)"
                )   << "mixture has only one component or number is not defined correctly, "
                    << "recheck numberOfComponents in the  thermophysicalProperties dict of your case"
                    << abort(FatalError);
            }

        label i;
        label j;
        label k=0;
 
        // CL: saves the number of mixtureCorrectionCoefficients needed for this mixture
        // CL: need to be set to 1
        label sizeOfVector_=1;

        //CL: size of the vector depends on the number of components
        //CL: determine the size of the vector    
        for (i=3; i<=nCom_;i++)
        {   
            sizeOfVector_=sizeOfVector_+(i-1);
        }
 
        //CL: setting the size 
        realMixtureCorrCoef_.setSize(sizeOfVector_);
       
        //CL: Reading the real mixture correction coefficients
        //CL: Naming convention e.g for 3 component mixture:
        //CL: mixtureCorrectionCoefficient_12, mixtureCorrectionCoefficient_13, mixtureCorrectionCoefficient_23
        for(i=1;i<=nCom_-1;i++)
        {
            for(j=i+1;j<=nCom_;j++)
            {
                realMixtureCorrCoef_[k]=dict.subDict("equationOfState")
                    .lookupOrDefault("mixtureCorrectionCoefficient_" + Foam::name(i) + Foam::name(j),0.0);
                k++;
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixtureSoaveRedlichKwong::write(Ostream& os) const
{
    soaveRedlichKwong::write(os);

}
*/
// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const mixtureSoaveRedlichKwong& srk)
{
    os  << static_cast<const soaveRedlichKwong&>(srk)<< token::SPACE ;

    os.check("Ostream& operator<<(Ostream& os, const mixtureSoaveRedlichKwong& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
