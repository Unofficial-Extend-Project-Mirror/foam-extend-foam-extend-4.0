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

Author
    Martin Beaudoin, Hydro-Quebec, 2014.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "controlSwitches.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace debug
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Type>
Foam::debug::controlSwitches<Type>::controlSwitches()
:
    switchValue_(Type(0))
{
}

template<class T>
Foam::debug::controlSwitches<T>::controlSwitches(const T& switchValue)
:
    switchValue_(switchValue)
{}


template<class T>
Foam::debug::controlSwitches<T>::controlSwitches
(
    const std::string& switchName,
    const T& switchValue,
    globalControlDictSwitchSet switchSet,
    std::map<std::string, std::list<controlSwitches<T> *> >** switchValuesTable
)
:
    switchSet_(switchSet),
    switchName_(switchName),
    switchValue_(switchValue)
{
    // Register the switch in its list
    if (*switchValuesTable == NULL)
    {
	*switchValuesTable = new std::map<std::string, std::list<controlSwitches<T> *> >();
    }

    switchValuesTable_ = *switchValuesTable;
    std::map<std::string, std::list<controlSwitches<T> *> >&switchValues = *switchValuesTable_;

    // Memorize this switch object address
    if (switchValues.find(switchName) != switchValues.end())
    {
	switchValues[switchName].push_back(this);
    }
    else
    {
	std::list<controlSwitches<T>* > pList;
	pList.push_back(this);
	switchValues.insert(std::pair<std::string, std::list<controlSwitches<T>* > >(switchName, pList));
    }
}


template<class T>
Foam::debug::controlSwitches<T>::controlSwitches(const Foam::debug::controlSwitches<T>& csw)
:
    switchName_(csw.switchName_),
    switchValue_(csw.switchValue_)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
Foam::debug::controlSwitches<T>::~controlSwitches()
{
    // Unregister the switch from its list
    std::map<std::string, std::list<controlSwitches<T> *> >&switchValuesTable = *switchValuesTable_;

    // Remove entry or key if pointers list is empty
    switchValuesTable[switchName_].remove(this);

    // Replace the updated list
    if(switchValuesTable[switchName_].size() == 0)
    {
	switchValuesTable.erase(switchName_);
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class T>
void Foam::debug::controlSwitches<T>::operator=(const Foam::debug::controlSwitches<T>& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
	std::cerr << "Foam::debug::controlSwitches<T>::operator=(const Foam::controlSwitches<T>&)"
		  << "--> FOAM FATAL ERROR: "
		  << "Attempted assignment to self"
		  << std::endl;
	std::abort();
	exit(-1);
    }
    else
    {
	switchValue_ = rhs.switchValue_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace debug

} // End namespace Foam

// ************************************************************************* //
