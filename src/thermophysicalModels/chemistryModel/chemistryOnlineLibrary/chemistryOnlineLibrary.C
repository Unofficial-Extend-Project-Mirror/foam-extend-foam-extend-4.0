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

#include "chemistryOnlineLibrary.H"
#include "error.H"
#include "chemistryModel.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::chemistryOnlineLibrary::staticData();


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construc from dictionary
Foam::chemistryOnlineLibrary::chemistryOnlineLibrary
(
    const chemistryModel& chemistry,
    const dictionary& chemistryProperties
)
:
    onlineDict_(chemistryProperties.subDict("onlineProperties")),
    chemistry_(chemistry),
    online_(onlineDict_.lookup("online")),
    chemisTree_(onlineDict_,*this),
    tolerance_(0.0),
    scaleFactor_(chemistry_.Y().size()+2),
    solutionScaleFactor_(chemistry_.Y().size()+2),
    logT_(false),
    clean_(false)
    
{

    if(online_)
    {
        tolerance_ = readScalar(onlineDict_.lookup("tolerance"));
        logT_ = onlineDict_.lookup("logT");
        clean_ = onlineDict_.lookup("cleanAll");
        dictionary scaleDict(onlineDict_.subDict("scaleFactor"));

        for(label i = 0; i<chemistry_.Y().size(); i++)
        {
            if(!scaleDict.found(chemistry_.Y()[i].name()))
            {
                scaleFactor_[i] = readScalar(scaleDict.lookup("otherSpecies"));
            }
            else
            {
                scaleFactor_[i] = readScalar(scaleDict.lookup(chemistry_.Y()[i].name()));
            }
        } 
        
        if(logT_)
        {
            scaleFactor_[chemistry_.Y().size()] = readScalar(scaleDict.lookup("logTemperature"));    
        }
        else
        {
            scaleFactor_[chemistry_.Y().size()] = readScalar(scaleDict.lookup("Temperature"));    
        }
        
        scaleFactor_[scaleFactor_.size() - 1] = readScalar(scaleDict.lookup("Pressure"));
        
        dictionary scaleSolutionDict(onlineDict_.subDict("scaleFactorSolution"));

        for(label i = 0; i<chemistry_.Y().size(); i++)
        {
            if(!scaleSolutionDict.found(chemistry_.Y()[i].name()))
            {
                solutionScaleFactor_[i] = readScalar(scaleSolutionDict.lookup("otherSpecies"));
            }
            else
            {
                solutionScaleFactor_[i] = readScalar(scaleSolutionDict.lookup(chemistry_.Y()[i].name()));
            }
        } 

        if(logT_)
        {
            solutionScaleFactor_[chemistry_.Y().size()] = readScalar(scaleSolutionDict.lookup("logTemperature"));    
        }
        else
        {
            solutionScaleFactor_[chemistry_.Y().size()] = readScalar(scaleSolutionDict.lookup("Temperature"));    
        }
        
        solutionScaleFactor_[solutionScaleFactor_.size() - 1] = readScalar(scaleSolutionDict.lookup("Pressure"));
        
               
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::chemistryOnlineLibrary> Foam::chemistryOnlineLibrary::New()
{
    return autoPtr<chemistryOnlineLibrary>(new chemistryOnlineLibrary);
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::chemistryOnlineLibrary::~chemistryOnlineLibrary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::chemPoint*  Foam::chemistryOnlineLibrary::findClosest
(
    const Foam::scalarField& v0
)
{
    return chemisTree_.findClosest(v0);
}


void Foam::chemistryOnlineLibrary::add
(
    const Foam::scalarField& v0,
    const Foam::scalarField& r,
    const Foam::scalarField& Wi,
    const Foam::scalar rhoi,
    const Foam::scalar& Tr,
    const Foam::scalar& pi,
    const Foam::scalar& deltaT
)
{

    scalarField vR = r;
    
    for(label j = 0; j<vR.size(); j++)
    {
        vR[j] *= Wi[j]/rhoi;
    }

/*
    vR.setSize(vR.size()+2);

    vR[vR.size()-2] = Tr;
    vR[vR.size()-1] = pi;

    chemPoint newPoint(v0, vR, scaleFactor_, solutionScaleFactor_, tolerance_, logT_, deltaT); 

*/
    
    vR.setSize(vR.size()+2);

    vR[vR.size()-2] = Tr;
    vR[vR.size()-1] = pi;
    
    scalarField vv0 = v0;
    
//    vv0.setSize(vv0.size()-1);
                   
    chemPoint newPoint(vv0, vR, scaleFactor_, solutionScaleFactor_, tolerance_, logT_, deltaT); 
    
    chemisTree_.insert(newPoint);    
    
}

const Foam::scalarField Foam::chemistryOnlineLibrary::calcNewC
(
    const Foam::chemPoint& c0
)
{
    
    return c0.r();    
    
}

void Foam::chemistryOnlineLibrary::clear()
{

    Info<< "Clearing chemistry library" << endl;

    //v_.clear();
    chemisTree_.cleanAll();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::chemistryOnlineLibrary::operator=(const Foam::chemistryOnlineLibrary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::chemistryOnlineLibrary::operator=(const Foam::chemistryOnlineLibrary&)")
            << "Attempted to assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
