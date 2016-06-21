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

Author
    Martin Beaudoin, Hydro-Quebec, 2014.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "controlSwitches.H"
#include <iomanip>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace debug
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::debug::controlSwitches<T>::controlSwitches
(
    const std::string& switchName,
    const T& switchValue,
    const std::string& switchDescription,
    globalControlDictSwitchSet switchSet,
    std::map<std::string, std::list<controlSwitches<T> *> >** switchValuesTable
)
:
    switchSet_(switchSet),
    switchName_(switchName),
    switchValue_(switchValue),
    switchDescription_(switchDescription)
{
    // Register the switch in its list
    if (*switchValuesTable == NULL)
    {
        *switchValuesTable =
            new std::map<std::string, std::list<controlSwitches<T> *> >();
    }

    switchValuesTable_ = *switchValuesTable;
    std::map<std::string, std::list<controlSwitches<T> *> >&
        switchValues = *switchValuesTable_;

    // Memorize this switch object address
    if (switchValues.find(switchName) != switchValues.end())
    {
        switchValues[switchName].push_back(this);
    }
    else
    {
        std::list<controlSwitches<T>* > pList;
        pList.push_back(this);

        switchValues.insert
        (
            std::pair<std::string, std::list<controlSwitches<T>* > >
            (
                switchName,
                pList
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
Foam::debug::controlSwitches<T>::~controlSwitches()
{
    // Unregister the switch from its list
    if (switchValuesTable_)
    {
        std::map<std::string, std::list<controlSwitches<T> *> >&
            switchValuesTable = *switchValuesTable_;

        // Remove entry or key if pointers list is empty
        switchValuesTable[switchName_].remove(this);

        // Replace the updated list
        if (switchValuesTable[switchName_].empty())
        {
            switchValuesTable.erase(switchName_);
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
void printControlSwitches
(
    const std::string& dictName,
    const std::map
    <
        std::string,
        std::list<Foam::debug::controlSwitches<T> *>
    >* mapListSwitchesPtr,
    const char* commsTypesNames[]
)
{
    std::cout << dictName << "\n{\n";

    if (mapListSwitchesPtr)
    {
        const std::map
        <
            std::string,
            std::list<Foam::debug::controlSwitches<T> *>
        >& mapListSwitches = *mapListSwitchesPtr;

        typename std::map
        <
            std::string,
            std::list<Foam::debug::controlSwitches<T> *>
        >::const_iterator it;

        // Compute the maximum length of the switches names.
        // Trying to make things pretty with std::setw by lining-up the values.
        // Useful for the debugSwitches because non-zero flags will be much
        // easier to spot.
        std::size_t maxLengthKey = 0;

        for (it = mapListSwitches.begin(); it != mapListSwitches.end(); it++)
        {
            std::string switchName = it->first;
            maxLengthKey = std::max(maxLengthKey, switchName.length());
        }

        // Still, we clip at 60 characters max.
        maxLengthKey = std::min(maxLengthKey, size_t(60));

        for (it = mapListSwitches.begin(); it != mapListSwitches.end(); it++)
        {
            // Switch name
            std::string switchName = it->first;

            // Switches list, all identical values, but from
            // different instances.
            // So we only need to echo the first switch of the list.
            std::list<Foam::debug::controlSwitches<T> *> swList = it->second;
            Foam::debug::controlSwitches<T>& value = *(swList.front());

            std::cout
                << "    " << std::left << std::setw(maxLengthKey)
                << switchName << " ";

            // Special handling for commsTypes from optimisationSwitches
            if (commsTypesNames && switchName.compare("commsType") == 0)
            {
                int valEnumAsIndex = value();

                std::cout << commsTypesNames[valEnumAsIndex] << ";";
            }
            else
            {
                std::cout << value() << ";";
            }

            // Now, for the switch description, since numerous switches might
            // be defined with identical names, but different descriptions
            // eg: ggi debugSwitch, we will concatenate all the non-empty
            // switches descriptions for a given switch

            std::string switchDescription("");

            typename
            std::list<Foam::debug::controlSwitches<T> *>::iterator itL;

            for (itL = swList.begin(); itL != swList.end(); itL++)
            {
                debug::controlSwitches<T>& sw = *(*itL);
                std::string thisSwitchDescr = sw.switchDescription();

                if
                (
                    !thisSwitchDescr.empty()
                 && switchDescription.find(thisSwitchDescr) ==
                    std::string::npos
                )
                {
                    switchDescription += thisSwitchDescr + ". ";
                }
            }

            if (!switchDescription.empty())
            {
                std::cout << "\t// " << switchDescription;
            }
            std::cout << std::endl;
        }
    }
    else
    {
        std::cout
            << "    // No switches of this type for this application"
            << std::endl;
    }

    std::cout << "}\n\n";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace debug


} // End namespace Foam

// ************************************************************************* //
