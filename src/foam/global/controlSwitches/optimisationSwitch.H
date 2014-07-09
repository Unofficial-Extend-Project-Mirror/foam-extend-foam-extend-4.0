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

Class
    Foam::infoSwitch

Author
    Martin Beaudoin, Hydro-Quebec, 2014.  All rights reserved

Description
    Handling of runtime controlSwitches.
    Specialization of controlSwitches class for optimisation switches

\*---------------------------------------------------------------------------*/

#ifndef optimisationSwitch_H
#define optimisationSwitch_H

#include "controlSwitches.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace debug
{

// Typedefs
// Hashtable to static class attributes addresses holding the
// runtime optimisation switch values
typedef std::map<std::string, std::list<controlSwitches<int> *> > ListOptimisationControlSwitches;
extern ListOptimisationControlSwitches* optimisationSwitchValues_;


class optimisationSwitch
    : public controlSwitches<int>
{
public:
    optimisationSwitch
    (
	const std::string& switchName,
	const int& switchValue
    )
	:
	controlSwitches<int>
	(
	    switchName,
	    debug::optimisationSwitchFromDict(switchName.c_str(), switchValue),
	    debug::OPTIMISATIONSWITCHES,
	    &optimisationSwitchValues_
	)
    {}


    // This constructor is essentially for 'commsType'
    optimisationSwitch
    (
	const std::string& switchName,
	const std::string& switchValue
    )
	:
	controlSwitches<int>
	(
	    switchName,
	    debug::optimisationSwitchFromDict(switchName, switchValue),
	    debug::OPTIMISATIONSWITCHES,
	    &optimisationSwitchValues_
	)
    {}

    virtual ~optimisationSwitch()
    {
    }

    // Assignement operator
    void operator=(const int& rhs)
    {
	controlSwitches<int>::operator=(rhs);
    }

    bool operator==(const int& rhs) const
    {
        return controlSwitches<int>::operator==(rhs);
    }

#if 0
    // Cannot use a conversion operator for Pstream::commsTypes
    // since including Pstream.H generates compilation errors.
    // MB 06/2014
    operator Pstream::commsTypes () const
    {
	return static_cast<Pstream::commsTypes>(operator()());
    }
#endif

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace debug


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
