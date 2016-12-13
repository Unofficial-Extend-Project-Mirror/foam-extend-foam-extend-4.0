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

Class
    Foam::debugSwitch

Author
    Martin Beaudoin, Hydro-Quebec, 2014.  All rights reserved

Description
    Handling of runtime controlSwitches.
    Specialization of controlSwitches class for debug switches

\*---------------------------------------------------------------------------*/

#ifndef debugSwitch_H
#define debugSwitch_H

#include "controlSwitches.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace debug
{

// Typedefs
// Hashtable to static class attributes addresses holding the
// runtime debug switch values
typedef std::map<std::string, std::list<controlSwitches<int> *> >
ListDebugControlSwitches;

extern ListDebugControlSwitches* debugSwitchValues_;

class debugSwitch
:
    public controlSwitches<int>
{
    // Private member functions

        //- Disallow construct as copy
        debugSwitch(const debugSwitch&);

        //- Disallow default bitwise assignment
        void operator=(const debugSwitch&);


public:

    // Constructors

        //- Construct from name and value
        debugSwitch
        (
            const std::string& switchName,
            const int& switchValue,
            const std::string& switchDescription = ""
        )
        :
            controlSwitches<int>
            (
                switchName,
                debug::debugSwitchFromDict(switchName.c_str(), switchValue),
                switchDescription,
                debug::DEBUG_SWITCHES,
                &debugSwitchValues_
            )
        {}


    //- Destructor
    virtual ~debugSwitch()
    {}


    // Member operators

        //- Assignement operator
        void operator=(const int& rhs)
        {
            controlSwitches<int>::operator=(rhs);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace debug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
