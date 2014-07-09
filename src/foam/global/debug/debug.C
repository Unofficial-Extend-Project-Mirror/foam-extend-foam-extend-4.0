/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Class for handling debugging switches.

\*---------------------------------------------------------------------------*/

#include "debug.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "dimensionedConstants.H"
#include "SortableList.H"
#include "debugSwitch.H"
#include "infoSwitch.H"
#include "optimisationSwitch.H"
#include "tolerancesSwitch.H"
#include "constantsSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace debug
{

//! @cond ignoreDocumentation - local scope
dictionary* controlDictPtr_(NULL);
dictionary* debugSwitchesPtr_(NULL);
dictionary* infoSwitchesPtr_(NULL);
dictionary* optimisationSwitchesPtr_(NULL);
dictionary* tolerancesPtr_(NULL);
dictionary* constantsPtr_(NULL);

// Hashtables to static class attributes addresses holding the
// runtime debug/info/optimisation/tolerances values
// This needs to go on the heap so the destructor will not get
//  called before the object's destructor it is overseeing
ListDebugControlSwitches* debugSwitchValues_(NULL);
ListInfoControlSwitches* infoSwitchValues_(NULL);
ListOptimisationControlSwitches* optimisationSwitchValues_(NULL);
ListTolerancesControlSwitches* tolerancesSwitchValues_(NULL);
ListConstantsControlSwitches* constantsSwitchValues_(NULL);


// to ensure controlDictPtr_ is deleted at the end of the run
class deleteControlDictPtr
{
public:

    deleteControlDictPtr()
    {}

    ~deleteControlDictPtr()
    {
        if (controlDictPtr_)
        {
            delete controlDictPtr_;
            controlDictPtr_ = 0;
        }
    }
};

deleteControlDictPtr deleteControlDictPtr_;
//! @endcond ignoreDocumentation


} // End namespace debug
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary& Foam::debug::controlDict()
{
    if (!controlDictPtr_)
    {
        // Allow users to override the location of the global controlDict
        // dictionary using an environment variable. Using this environment
        // variable, one can assign a different global controlDict for each
        // case, without having to modify the "default" ones.
        fileName globControlDictFileName = getEnv("FOAM_GLOBAL_CONTROLDICT");

        // Fallback to default locations if filename is empty or not valid
        if( ! isFile(globControlDictFileName) )
            globControlDictFileName = findEtcFile("controlDict", true);

        controlDictPtr_ = new dictionary
        (
            IFstream(globControlDictFileName)()
        );
    }

    return *controlDictPtr_;
}


Foam::dictionary& Foam::debug::switchSet
(
    const char* subDictName,
    dictionary*& subDictPtr
)
{
    if (!subDictPtr)
    {
        entry* ePtr = controlDict().lookupEntryPtr
        (
            subDictName, false, false
        );

        if (!ePtr || !ePtr->isDict())
        {
            cerr<< "debug::switchSet(const char*, dictionary*&):\n"
                << "    Cannot find " <<  subDictName << " in dictionary "
                << controlDict().name().c_str()
                << std::endl << std::endl;

            ::exit(1);
        }

        subDictPtr = &ePtr->dict();
    }

    return *subDictPtr;
}


Foam::dictionary& Foam::debug::debugSwitches()
{
    return switchSet("DebugSwitches", debugSwitchesPtr_);
}


Foam::dictionary& Foam::debug::infoSwitches()
{
    return switchSet("InfoSwitches", infoSwitchesPtr_);
}


Foam::dictionary& Foam::debug::optimisationSwitches()
{
    return switchSet("OptimisationSwitches", optimisationSwitchesPtr_);
}


Foam::dictionary& Foam::debug::tolerances()
{
    return switchSet("Tolerances", tolerancesPtr_);
}

Foam::dictionary& Foam::debug::constants()
{
    return switchSet("DimensionedConstants", constantsPtr_);
}


int Foam::debug::debugSwitchFromDict(const char* name, const int defaultValue)
{
    return debugSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


int Foam::debug::infoSwitchFromDict(const char* name, const int defaultValue)
{
    return infoSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


int Foam::debug::optimisationSwitchFromDict
(
    const char* name,
    const int defaultValue
)
{
    return optimisationSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


// This is essentially for reading 'commsType'
int Foam::debug::optimisationSwitchFromDict
(
    const Foam::string name,
    const Foam::word defaultValueStr
)
{
    int retValue = 0;

    if (name == "commsType")
    {
	Foam::string valueStr = optimisationSwitches().lookupOrAddDefault
	(
	    name, defaultValueStr, false, false
	);

	retValue = Pstream::commsTypeNames[valueStr];
    }
    else
    {
	// Need a warning here...
    }

    return retValue;
}


double Foam::debug::tolerancesFromDict
(
    const char* name,
    const double defaultValue
)
{
    return tolerances().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}

double Foam::debug::constantsFromDict
(
    const char* name,
    const double defaultValue
)
{
    return tolerances().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}

void Foam::debug::updateCentralDictVars
(
    Foam::debug::globalControlDictSwitchSet globalControlDictSwitchSetName,
    Foam::string keyValues
)
{
    Foam::word token;

    // Extract the multiple keyValues, separated by ','
    std::istringstream issmkvs(keyValues);
    while ( getline(issmkvs, token, ',') )
    {
	// Cleanup
	Foam::string::stripInvalid<word>(token);

	// Extract the keyValue pair, separated by '='
	std::istringstream isskv(token);

	while ( getline(isskv, token, '=') )
	{
	    Foam::word key = token;

	    switch(globalControlDictSwitchSetName)
	    {
	        case DEBUGSWITCHES:
		{
		    if (debugSwitchValues_)
		    {
			ListDebugControlSwitches& debugSwitchValues =
			    *debugSwitchValues_;

			label oldDebugValue;
			label newDebugValue;
			isskv >> newDebugValue;

			if (Foam::debug::debugSwitches().readIfPresent
			    (
				key,
				oldDebugValue
			    )
			   )
			{
			    Info << endl
				 << "Warning: Modification of DebugSwitch: "
				 << key << endl
				 << "    Old value: " << oldDebugValue << endl
				 << "    New value: " << newDebugValue << endl
				 << endl;

			    Foam::debug::debugSwitches().set
				(
				    key,
				    newDebugValue
				);

			    std::list<controlSwitches<int> *> curList =
				debugSwitchValues[key];

			    // Modify all entries for this key
			    forAllIter
			    (
				std::list<controlSwitches<int> *>,
				curList,
				iterI
			    )
			    {
				*(*iterI) =  newDebugValue;
			    }
			}
			else
			{
			    //  Usage of non-existent DebugSwitches: best to
			    // abort right away
			    SortableList<Foam::word> sortedValidKeys
			    (
				debugSwitchValues.size()
			    );

			    int i=0;
			    forAllIter
			    (
				ListDebugControlSwitches,
				debugSwitchValues,
				iterI
			    )
			    {
				sortedValidKeys[i++] = iterI->first;
			    }
			    sortedValidKeys.sort();

			    FatalError
				<< "Usage of non-existent DebugSwitches name: "
				<< key
				<< endl << endl
				<< "Valid entries for this application are: "
				<< sortedValidKeys
				<< exit(FatalError);
			}
		    }
		    else
		    {
			FatalError
			    << "No DebugSwitches values are available for "
			    << "this application."
			    << exit(FatalError);
		    }
		}
		break;
	        case INFOSWITCHES:
		{
		    if (infoSwitchValues_)
		    {
			ListInfoControlSwitches& infoSwitchValues =
			    *infoSwitchValues_;

			label oldInfoValue;
			label newInfoValue;
			isskv >> newInfoValue;

			if (Foam::debug::infoSwitches().readIfPresent
			    (
				key,
				oldInfoValue
			    )
			   )
			{
			    Info << endl
				 << "Warning: Modification of InfoSwitch: "
				 << key << endl
				 << "    Old value: " << oldInfoValue << endl
				 << "    New value: " << newInfoValue << endl
				 << endl;

			    Foam::debug::infoSwitches().set(key, newInfoValue);

			    std::list<controlSwitches<int> *> curList =
				infoSwitchValues[key];

			    // Modify all entries for this key
			    forAllIter
			    (
				std::list<controlSwitches<int> *>,
				curList,
				iterI
			    )
			    {
				*(*iterI) =  newInfoValue;
			    }
			}
			else
			{
			    // Usage of non-existent InfoSwitches: best to
			    // abort right away
			    SortableList<Foam::word> sortedValidKeys
			    (
				infoSwitchValues.size()
			    );

			    int i=0;
			    forAllIter
			    (
				ListInfoControlSwitches,
				infoSwitchValues,
				iterI
			    )
			    {
				sortedValidKeys[i++] = iterI->first;
			    }
			    sortedValidKeys.sort();

			    FatalError
				<< "Usage of non-existent InfoSwitches name: "
				<< key
				<< endl << endl
				<< "Valid entries for this application are: "
				<< sortedValidKeys
				<< exit(FatalError);
			}
		    }
		    else
		    {
			FatalError
			    << "No InfoSwitches values are available for "
			    << "this application."
			    << exit(FatalError);
		    }
		}
		break;
	        case OPTIMISATIONSWITCHES:
		{
		    if (optimisationSwitchValues_)
		    {
			ListOptimisationControlSwitches&
			    optimisationSwitchValues = *optimisationSwitchValues_;

			label newOptimisationValue;
			label oldOptimisationValue;
			word newOptimisationValueStr;
			word oldOptimisationValueStr;
			bool keyIsPresent(false);

			// We need to check if the switch 'commsType' is being
			// overriden. This switch is coded as an enum, but will
			// be specified using a string like "blocking",
			// "scheduled, etc. Some additional logic is needed.
			if (key == "commsType")
			{
			    // Handle a string value, then convert to enum
			    isskv >> newOptimisationValueStr;
			    newOptimisationValue =
				Pstream::commsTypeNames[newOptimisationValueStr];
			    keyIsPresent =
				Foam::debug::optimisationSwitches().readIfPresent
				(
				    key,
				    oldOptimisationValueStr
				);
			}
			else
			{
			    // Handling label values
			    isskv >> newOptimisationValue;
			    keyIsPresent =
				Foam::debug::optimisationSwitches().readIfPresent
				(
				    key,
				    oldOptimisationValue
				);

			    std::ostringstream newOptimisationValueOstream;
			    newOptimisationValueOstream << newOptimisationValue;
			    newOptimisationValueStr =
				newOptimisationValueOstream.str();

			    std::ostringstream oldOptimisationValueOstream;
			    oldOptimisationValueOstream << oldOptimisationValue;
			    oldOptimisationValueStr =
				oldOptimisationValueOstream.str();
			}

			if(keyIsPresent)
			{
			    Info << endl
				 << "Warning: Modification of "
				 << "OptimisationSwitch: "
				 << key << endl
				 << "    Old value: "
				 << oldOptimisationValueStr << endl
				 << "    New value: "
				 << newOptimisationValueStr << endl
				 << endl;

			    Foam::debug::optimisationSwitches().set
			    (
				key,
				newOptimisationValue
			    );

			    std::list<controlSwitches<int> *> curList =
				optimisationSwitchValues[key];

			    std::cout << "curList.size(): "  << curList.size() << std::endl;
			    // Modify all entries for this key
			    forAllIter
			    (
				std::list<controlSwitches<int> *>,
				curList,
				iterI
			    )
			    {
				*(*iterI) =  newOptimisationValue;
			    }
			}
			else
		        {
			    // Usage of non-existent OptimisationSwitches: best to
			    // abort right away
			    SortableList<Foam::word> sortedValidKeys
			    (
				optimisationSwitchValues.size()
			    );

			    int i=0;
			    forAllIter
			    (
				ListOptimisationControlSwitches,
				optimisationSwitchValues,
				iterI
			    )
			    {
				sortedValidKeys[i++] = iterI->first;
			    }
			    sortedValidKeys.sort();

			    FatalError
				<< "Usage of non-existent "
				<< " OptimisationSwitches name: " << key
				<< endl << endl
				<< "Valid entries for this application are: "
				<< sortedValidKeys
				<< exit(FatalError);
			}
		    }
		    else
		    {
			FatalError
			    << "No OptimisationSwitches values are available "
			    << "for this application."
			    << exit(FatalError);
		    }
		}
		break;
	        case TOLERANCES:
		{
		    if (tolerancesSwitchValues_)
		    {
			ListTolerancesControlSwitches& tolerancesSwitchValues =
			    *tolerancesSwitchValues_;

			scalar oldTolerancesValue;
			scalar newTolerancesValue;
			isskv >> newTolerancesValue;

			if (Foam::debug::tolerances().readIfPresent
			    (
				key,
				oldTolerancesValue
			    )
			   )
			{
			    Info << endl
				 << "Warning: Modification of Tolerance: "
				 << key << endl
				 << "    Old value: "
				 << oldTolerancesValue << endl
				 << "    New value: "
				 << newTolerancesValue << endl
				 << endl;

			    Foam::debug::tolerances().set
			    (
				key,
				newTolerancesValue
			    );

			    std::list<controlSwitches<scalar> *> curList =
				tolerancesSwitchValues[key];

			    // Modify all entries for this key
			    forAllIter
			    (
				std::list<controlSwitches<scalar> *>,
				curList,
				iterI
			    )
			    {
				*(*iterI) =  newTolerancesValue;
			    }
			}
			else
			{
			    // Usage of non-existent Tolerances: best to
			    // abort right away
			    SortableList<Foam::word> sortedValidKeys
			    (
				tolerancesSwitchValues.size()
			    );

			    int i=0;
			    forAllIter
			    (
				ListTolerancesControlSwitches,
				tolerancesSwitchValues,
				iterI
			    )
			    {
				sortedValidKeys[i++] = iterI->first;
			    }
			    sortedValidKeys.sort();

			    FatalError
				<< "Usage of non-existent Tolerances name: "
				<< key
				<< endl << endl
				<< "Valid entries for this application are: "
				<< sortedValidKeys
				<< exit(FatalError);
			}
		    }
		    else
		    {
			FatalError
			    << "No Tolerances values are available for "
			    << "this application."
			    << exit(FatalError);
		    }
		}
		break;
	        case DIMENSIONEDCONSTANTS:
		{
		    if (constantsSwitchValues_)
		    {
			ListConstantsControlSwitches& constantsSwitchValues =
			    *constantsSwitchValues_;

			scalar oldDimensionedConstantsValue;
			scalar newDimensionedConstantsValue;
			isskv >> newDimensionedConstantsValue;

			if (Foam::dimensionedConstants().readIfPresent
			    (
				key,
				oldDimensionedConstantsValue
			    )
			  )
			{
			    Info << endl
				 << "Warning: Modification of DimensionedConstant: "
				 << key << endl
				 << "    Old value: "
				 << oldDimensionedConstantsValue << endl
				 << "    New value: "
				 << newDimensionedConstantsValue << endl
				 << endl;

			    Foam::dimensionedConstants().set
			    (
				key,
				newDimensionedConstantsValue
			    );

			    std::list<controlSwitches<scalar> *> curList =
				constantsSwitchValues[key];

			    // Modify all entries for this key
			    forAllIter
			    (
				std::list<controlSwitches<scalar> *>,
				curList,
				iterI
			    )
			    {
				*(*iterI) =  newDimensionedConstantsValue;
			    }
			}
			else
			{
			    // Usage of non-existent DimensionedConstants: best to
			    // abort right away
			    SortableList<Foam::word> sortedValidKeys
			    (
				constantsSwitchValues.size()
			    );

			    int i=0;
			    forAllIter
			    (
				ListConstantsControlSwitches,
				constantsSwitchValues,
				iterI
			    )
			    {
				sortedValidKeys[i++] = iterI->first;
			    }
			    sortedValidKeys.sort();

			    FatalError
				<< "Usage of non-existent "
				<< "DimensionedConstants name: "
				<< key
				<< endl << endl
				<< "Valid entries for this application are: "
				<< sortedValidKeys
				<< exit(FatalError);
			}
		    }
		    else
		    {
			FatalError
			    << "No DimensionedConstants values are available "
			    << "for this application."
			    << exit(FatalError);
		    }
		}
		break;
	        default:
		break;
	    }
	}
    }
}

// ************************************************************************* //
