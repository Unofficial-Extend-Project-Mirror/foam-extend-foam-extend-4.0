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
    Class for handling debugging switches.

\*---------------------------------------------------------------------------*/

#include "debug.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace debug
{

dictionary* controlDictPtr_(NULL);
dictionary* debugSwitchesPtr_(NULL);
dictionary* infoSwitchesPtr_(NULL);
dictionary* optimisationSwitchesPtr_(NULL);
dictionary* tolerancesPtr_(NULL);

//- Class to ensure controlDictPtr_ is deleted at the end of the run
//  @cond ignore documentation for this class
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
        }
    }
};
//! @endcond

deleteControlDictPtr deleteControlDictPtr_;


dictionary& switchSet(const char* switchSetName, dictionary* switchSetPtr)
{
    if (!switchSetPtr)
    {
        if (!controlDict().found(switchSetName))
        {
            cerr<< "debug::switchSet(const char*, dictionary*): " << std::endl
                << "    Cannot find " <<  switchSetName
            << " in dictionary " << controlDictPtr_->name().c_str()
            << std::endl << std::endl;

            ::exit(1);
        }

        switchSetPtr =
            const_cast<dictionary*>(&(controlDict().subDict(switchSetName)));
    }

    return *switchSetPtr;
}


int debugSwitch
(
    dictionary& switchSet,
    const char* switchName,
    const int defaultValue
)
{
    if (switchSet.found(switchName))
    {
        return readInt(switchSet.lookup(switchName));
    }
    else
    {
        switchSet.add(switchName, defaultValue);
        return defaultValue;
    }
}


double debugTolerance
(
    dictionary& switchSet,
    const char* switchName,
    const double defaultValue
)
{
    if (switchSet.found(switchName))
    {
        return readDoubleScalar(switchSet.lookup(switchName));
    }
    else
    {
        switchSet.add(switchName, defaultValue);
        return defaultValue;
    }
}

}  // End namespace debug


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary& debug::controlDict()
{
    if (!controlDictPtr_)
    {
        fileName controlDictFileName(dotFoam("controlDict"));

        IFstream dictFile(controlDictFileName);

        if (!dictFile.good())
        {
            cerr<< "debug::controlDict(): "
                << "Cannot open essential file " << controlDictFileName.c_str()
                << std::endl << std::endl;
            ::exit(1);
        }

        controlDictPtr_ = new dictionary(dictFile);
    }

    return *controlDictPtr_;
}


dictionary& debug::debugSwitches()
{
    return switchSet("DebugSwitches", debugSwitchesPtr_);
}


int debug::debugSwitch(const char* switchName, const int defaultValue)
{
    return debugSwitch
    (
        debugSwitches(),
        switchName,
        defaultValue
    );
}


dictionary& debug::infoSwitches()
{
    return switchSet("InfoSwitches", infoSwitchesPtr_);
}


int debug::infoSwitch(const char* switchName, const int defaultValue)
{
    return debugSwitch
    (
        infoSwitches(),
        switchName,
        defaultValue
    );
}


dictionary& debug::optimisationSwitches()
{
    return switchSet("OptimisationSwitches", optimisationSwitchesPtr_);
}


int debug::optimisationSwitch(const char* switchName, const int defaultValue)
{
    return debugSwitch
    (
        optimisationSwitches(),
        switchName,
        defaultValue
    );
}


dictionary& debug::tolerances()
{
    return switchSet("Tolerances", tolerancesPtr_);
}


double debug::tolerances(const char* switchName, const double defaultValue)
{
    return debugTolerance
    (
        tolerances(),
        switchName,
        defaultValue
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
