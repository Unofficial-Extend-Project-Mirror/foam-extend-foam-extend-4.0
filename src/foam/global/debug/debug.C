/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
#include "fileName.H"
#include "NamedEnum.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char*
Foam::NamedEnum
<
    Foam::debug::globalControlDictSwitchSet,
    Foam::debug::DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
>::names[] =
{
    "DebugSwitches",
    "InfoSwitches",
    "OptimisationSwitches",
    "Tolerances",
    "DimensionedConstants"
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace debug
{

//! @cond ignoreDocumentation - local scope
dictionary* controlDictPtr_(NULL);

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
            controlDictPtr_ = NULL;
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
        // We are totally getting rid of the central controlDict file located
        // under $WM_PROJECT_DIR/etc/controlDict.
        // But we still need to create a central dictionary for
        // global control switches.
        // At first, all the different subdicts will be empty, like this:
        //     DebugSwitches{}
        //     InfoSwitches{}
        //     OptimisationSwitches{}
        //     Tolerances{}
        //     DimensionedConstants{}
        //
        // We will also merge with other predefined dictionaries
        // as we go along.
        // Then, as the libraries will instantiate the control switches,
        // the corresponding subdicts will get populated.
        // MB 05/2015

        // Make control dict

        controlDictPtr_ = new dictionary();
        dictionary& cDict = *controlDictPtr_;

        // Cannot use static NamedEnum in static initialisation
        // Making local copy
        const NamedEnum
        <
            globalControlDictSwitchSet,
            DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
        > globalControlDictSwitchSetNames;

        // Add a subDict for each switch set
        forAll (globalControlDictSwitchSetNames, gI)
        {
            cDict.add(globalControlDictSwitchSetNames.names[gI], dictionary());
        }
    }

    return *controlDictPtr_;
}


Foam::dictionary& Foam::debug::debugSwitches()
{
    // Cannot use static NamedEnum in static initialisation
    // Making local copy
    const NamedEnum
    <
        globalControlDictSwitchSet,
        DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
    > globalControlDictSwitchSetNames;

    return controlDict().subDict
    (
        globalControlDictSwitchSetNames[debug::DEBUG_SWITCHES]
    );
}


Foam::dictionary& Foam::debug::infoSwitches()
{
    // Cannot use static NamedEnum in static initialisation
    // Making local copy
    const NamedEnum
    <
        globalControlDictSwitchSet,
        DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
    > globalControlDictSwitchSetNames;

    return controlDict().subDict
    (
        globalControlDictSwitchSetNames[debug::INFO_SWITCHES]
    );
}


Foam::dictionary& Foam::debug::optimisationSwitches()
{
    // Cannot use static NamedEnum in static initialisation
    // Making local copy
    const NamedEnum
    <
        globalControlDictSwitchSet,
        DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
    > globalControlDictSwitchSetNames;

    return controlDict().subDict
    (
        globalControlDictSwitchSetNames[debug::OPTIMISATION_SWITCHES]
    );
}


Foam::dictionary& Foam::debug::tolerances()
{
    // Cannot use static NamedEnum in static initialisation
    // Making local copy
    const NamedEnum
    <
        globalControlDictSwitchSet,
        DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
    > globalControlDictSwitchSetNames;

    return controlDict().subDict
    (
        globalControlDictSwitchSetNames[debug::TOLERANCES]
    );
}


Foam::dictionary& Foam::debug::constants()
{
    // Cannot use static NamedEnum in static initialisation
    // Making local copy
    const NamedEnum
    <
        globalControlDictSwitchSet,
        DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
    > globalControlDictSwitchSetNames;

    return controlDict().subDict
    (
        globalControlDictSwitchSetNames[debug::DIMENSIONED_CONSTANTS]
    );
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
    const string name,
    const word defaultValueStr
)
{
    int retValue = 0;

    if (name == "commsType")
    {
        string valueStr = optimisationSwitches().lookupOrAddDefault
        (
            name,
            defaultValueStr,
            false,
            false
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
    return constants().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}

void Foam::debug::updateCentralDictVars
(
    const debug::globalControlDictSwitchSet gcdSwitchSetName,
    const string& keyValues,
    const bool verbose
)
{
    word token;

    // Extract the multiple keyValues, separated by ','
    std::istringstream issmkvs(keyValues);

    while (getline(issmkvs, token, ','))
    {
        // Cleanup
        string::stripInvalid<word>(token);

        // Extract the keyValue pair, separated by '='
        std::istringstream isskv(token);

        while (getline(isskv, token, '='))
        {
            word key = token;

            switch(gcdSwitchSetName)
            {
                case DEBUG_SWITCHES:
                {
                    if (debugSwitchValues_)
                    {
                        ListDebugControlSwitches& debugSwitchValues =
                            *debugSwitchValues_;

                        label oldDebugValue;
                        label newDebugValue;
                        isskv >> newDebugValue;

                        if
                        (
                            debug::debugSwitches().readIfPresent
                            (
                                key,
                                oldDebugValue
                            )
                        )
                        {
                            if (verbose)
                            {
                                InfoIn
                                (
                                    "debug::updateCentralDictVars\n"
                                    "(\n"
                                    "    const globalControlDictSwitchSet,\n"
                                    "    const string& keyValues,\n"
                                    "    const bool verbose\n"
                                    ")"
                                )   << "CtrlDict : "
                                    << "Modification of DebugSwitches: "
                                    << key
                                    << ": old value: " << oldDebugValue
                                    << ", new value: " << newDebugValue
                                    << endl;
                            }

                            debug::debugSwitches().set
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
                            SortableList<word> sortedValidKeys
                            (
                                debugSwitchValues.size()
                            );

                            int i = 0;
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

                            WarningIn
                            (
                                "debug::updateCentralDictVars\n"
                                "(\n"
                                "    const globalControlDictSwitchSet,\n"
                                "    const string& keyValues,\n"
                                "    const bool verbose\n"
                                ")"
                            )   << "Usage of non-existent DebugSwitches name: "
                                << key << ".  Cannot be changed." << nl << nl
                                << "Valid entries for this application are: "
                                << sortedValidKeys
                                << endl;
                        }
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "debug::updateCentralDictVars\n"
                            "(\n"
                            "    const globalControlDictSwitchSet,\n"
                            "    const string& keyValues,\n"
                            "    const bool verbose\n"
                            ")"
                        )   << "No DebugSwitches values are available for "
                            << "this application."
                            << exit(FatalError);
                    }
                }
                break;
                case INFO_SWITCHES:
                {
                    if (infoSwitchValues_)
                    {
                        ListInfoControlSwitches& infoSwitchValues =
                            *infoSwitchValues_;

                        label oldInfoValue;
                        label newInfoValue;
                        isskv >> newInfoValue;

                        if
                        (
                            debug::infoSwitches().readIfPresent
                            (
                                key,
                                oldInfoValue
                            )
                        )
                        {
                            if (verbose)
                            {
                                Info<< "CtrlDict : "
                                    << "Modification of InfoSwitches: "
                                    << key
                                    << ": old value: " << oldInfoValue
                                    << ", new value: " << newInfoValue
                                    << endl;
                            }

                            debug::infoSwitches().set(key, newInfoValue);

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
                            SortableList<word> sortedValidKeys
                            (
                                infoSwitchValues.size()
                            );

                            int i = 0;
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

                            WarningIn
                            (
                                "debug::updateCentralDictVars\n"
                                "(\n"
                                "    const globalControlDictSwitchSet,\n"
                                "    const string& keyValues,\n"
                                "    const bool verbose\n"
                                ")"
                            )   << "Usage of non-existent InfoSwitches name: "
                                << key << ".  Cannot change." << nl << nl
                                << "Valid entries for this application are: "
                                << sortedValidKeys
                                << endl;
                        }
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "debug::updateCentralDictVars\n"
                            "(\n"
                            "    const globalControlDictSwitchSet,\n"
                            "    const string& keyValues,\n"
                            "    const bool verbose\n"
                            ")"
                        )   << "No InfoSwitches values are available for "
                            << "this application."
                            << exit(FatalError);
                    }
                }
                break;
                case OPTIMISATION_SWITCHES:
                {
                    if (optimisationSwitchValues_)
                    {
                        ListOptimisationControlSwitches&
                            optimisationSwitchValues =
                            *optimisationSwitchValues_;

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
                                Pstream::commsTypeNames
                                [
                                    newOptimisationValueStr
                                ];

                            keyIsPresent =
                                debug::optimisationSwitches().readIfPresent
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
                                debug::optimisationSwitches().readIfPresent
                                (
                                    key,
                                    oldOptimisationValue
                                );

                            std::ostringstream newOptimisationValueOstream;

                            newOptimisationValueOstream
                                << newOptimisationValue;

                            newOptimisationValueStr =
                                newOptimisationValueOstream.str();

                            std::ostringstream oldOptimisationValueOstream;

                            oldOptimisationValueOstream
                                << oldOptimisationValue;

                            oldOptimisationValueStr =
                                oldOptimisationValueOstream.str();
                        }

                        if (keyIsPresent)
                        {
                            if (verbose)
                            {
                                InfoIn
                                (
                                    "\ndebug::updateCentralDictVars\n"
                                    "(\n"
                                    "    const globalControlDictSwitchSet,\n"
                                    "    const string& keyValues,\n"
                                    "    const bool verbose\n"
                                    ")"
                                )   << "CtrlDict : Modification of "
                                    << "OptimisationSwitches: "
                                    << key
                                    << ": old value: "
                                    << oldOptimisationValueStr
                                    << ", new value: "
                                    << newOptimisationValueStr
                                    << endl;
                            }

                            debug::optimisationSwitches().set
                            (
                                key,
                                newOptimisationValue
                            );

                            std::list<controlSwitches<int> *> curList =
                                optimisationSwitchValues[key];

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
                            // Usage of non-existent OptimisationSwitches:
                            // best to abort right away
                            SortableList<word> sortedValidKeys
                            (
                                optimisationSwitchValues.size()
                            );

                            int i = 0;
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

                            WarningIn
                            (
                                "debug::updateCentralDictVars\n"
                                "(\n"
                                "    const globalControlDictSwitchSet,\n"
                                "    const string& keyValues,\n"
                                "    const bool verbose\n"
                                ")"
                            )   << "Usage of non-existent "
                                << " OptimisationSwitches name: " << key
                                << ".  Cannot change." << endl << endl
                                << "Valid entries for this application are: "
                                << sortedValidKeys
                                << endl;
                        }
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "debug::updateCentralDictVars\n"
                            "(\n"
                            "    const globalControlDictSwitchSet,\n"
                            "    const string& keyValues,\n"
                            "    const bool verbose\n"
                            ")"
                        )   << "No OptimisationSwitches values are available "
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

                        if
                        (
                            debug::tolerances().readIfPresent
                            (
                                key,
                                oldTolerancesValue
                            )
                        )
                        {
                            if (verbose)
                            {
                                InfoIn
                                (
                                    "debug::updateCentralDictVars\n"
                                    "(\n"
                                    "    const globalControlDictSwitchSet,\n"
                                    "    const string& keyValues,\n"
                                    "    const bool verbose\n"
                                    ")"
                                )   << "CtrlDict : "
                                    << "Modification of Tolerances: "
                                    << key
                                    << ": old value: " << oldTolerancesValue
                                    << ", new value: " << newTolerancesValue
                                    << endl;
                            }

                            debug::tolerances().set
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
                            SortableList<word> sortedValidKeys
                            (
                                tolerancesSwitchValues.size()
                            );

                            int i = 0;
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

                            WarningIn
                            (
                                "debug::updateCentralDictVars\n"
                                "(\n"
                                "    const globalControlDictSwitchSet,\n"
                                "    const string& keyValues,\n"
                                "    const bool verbose\n"
                                ")"
                            )   << "Usage of non-existent Tolerances name: "
                                << key << ".  Cannot change." << nl << nl
                                << "Valid entries for this application are: "
                                << sortedValidKeys
                                << endl;
                        }
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "debug::updateCentralDictVars\n"
                            "(\n"
                            "    const globalControlDictSwitchSet,\n"
                            "    const string& keyValues,\n"
                            "    const bool verbose\n"
                            ")"
                        )   << "No Tolerances values are available for "
                            << "this application."
                            << exit(FatalError);
                    }
                }
                break;
                case DIMENSIONED_CONSTANTS:
                {
                    if (constantsSwitchValues_)
                    {
                        ListConstantsControlSwitches& constantsSwitchValues =
                            *constantsSwitchValues_;

                        scalar oldDimensionedConstantsValue;
                        scalar newDimensionedConstantsValue;
                        isskv >> newDimensionedConstantsValue;

                        if
                        (
                            dimensionedConstants().readIfPresent
                            (
                                key,
                                oldDimensionedConstantsValue
                            )
                        )
                        {
                            if (verbose)
                            {
                                InfoIn
                                (
                                    "debug::updateCentralDictVars\n"
                                    "(\n"
                                    "    const globalControlDictSwitchSet,\n"
                                    "    const string& keyValues,\n"
                                    "    const bool verbose\n"
                                    ")"
                                )   << "CtrlDict : "
                                    << "Modification of DimensionedConstants: "
                                    << key
                                    << ": old value: "
                                    << oldDimensionedConstantsValue
                                    << ", new value: "
                                    << newDimensionedConstantsValue
                                    << endl;
                            }

                            dimensionedConstants().set
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
                            // Usage of non-existent DimensionedConstants:
                            // best to abort right away
                            SortableList<word> sortedValidKeys
                            (
                                constantsSwitchValues.size()
                            );

                            int i = 0;
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

                            WarningIn
                            (
                                "debug::updateCentralDictVars\n"
                                "(\n"
                                "    const globalControlDictSwitchSet,\n"
                                "    const string& keyValues,\n"
                                "    const bool verbose\n"
                                ")"
                            )   << "Usage of non-existent "
                                << "DimensionedConstants name: "
                                << key << "Cannot change." << nl << nl
                                << "Valid entries for this application are: "
                                << sortedValidKeys
                                << endl;
                        }
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "debug::updateCentralDictVars\n"
                            "(\n"
                            "    const globalControlDictSwitchSet,\n"
                            "    const string& keyValues,\n"
                            "    const bool verbose\n"
                            ")"
                        )   << "No DimensionedConstants values are available "
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


void Foam::debug::updateCentralDictVars
(
    const fileName& fName,
    const bool verbose
)
{
    // Add a trace to the banner
    if (verbose)
    {
        Info << "CtrlDict : " << fName << endl;
    }

    if (fName.size() && isFile(fName))
    {
        const dictionary* optionalGlobControlDictPtr_ = new dictionary
        (
            IFstream(fName)()
        );

        const dictionary& optionalGlobControlDict =
            *optionalGlobControlDictPtr_;

        const NamedEnum
        <
            debug::globalControlDictSwitchSet,
            DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
        >
        gcdSwitchSetNames;

        // We merge the content of each subdicts
        forAll(gcdSwitchSetNames, gI)
        {
            string subDictName(gcdSwitchSetNames.names[gI]);

            if (optionalGlobControlDict.isDict(subDictName))
            {
                const dictionary& subD =
                    optionalGlobControlDict.subDict(subDictName);
                const wordList toc = subD.toc();

                OStringStream keyValuePairs;
                forAll(toc, tocI)
                {
                    token value;
                    subD.lookup(toc[tocI]).read(value);
                    keyValuePairs << toc[tocI] << "=" << value << ", ";
                }

                // Update global dictionary
                updateCentralDictVars
                (
                    gcdSwitchSetNames[subDictName],
                    keyValuePairs.str(),
                    false // keep quiet when using files...
                );
            }
        }

        // Cleaning up
        delete optionalGlobControlDictPtr_;
    }
}

void Foam::debug::dumpControlSwitchesToConsole()
{
    const NamedEnum
    <
        debug::globalControlDictSwitchSet,
        DIM_GLOBAL_CONTROL_DICT_SWITCH_SET
    >
    globalControlDictSwitchSetNames;

    Info << endl;

    debug::printControlSwitches
    (
        globalControlDictSwitchSetNames[debug::DEBUG_SWITCHES],
        debug::debugSwitchValues_
    );

    debug::printControlSwitches
    (
        globalControlDictSwitchSetNames[debug::INFO_SWITCHES],
        debug::infoSwitchValues_
    );

    // We are forced to pass the string descriptions of the
    // Pstream::commsTypes for the optimisationSwitches group because
    // this switch is in fact an enum but we need to specify its
    // corresponding string equivalent in a controlDict
    // dictionary. And at the low level we are playing, including
    // Pstream.H is out of the question.  MB 2015
    debug::printControlSwitches
    (
        globalControlDictSwitchSetNames[debug::OPTIMISATION_SWITCHES],
        debug::optimisationSwitchValues_,
        Pstream::commsTypeNames.names
    );

    debug::printControlSwitches
    (
        globalControlDictSwitchSetNames[debug::TOLERANCES],
        debug::tolerancesSwitchValues_
    );

    debug::printControlSwitches
    (
        globalControlDictSwitchSetNames[debug::DIMENSIONED_CONSTANTS],
        debug::constantsSwitchValues_
    );

    Info << endl;

    return;
}


// ************************************************************************* //
