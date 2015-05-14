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
    Specialization of controlSwitches class for info switches

\*---------------------------------------------------------------------------*/

#ifndef infoSwitch_H
#define infoSwitch_H

#include "controlSwitches.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace debug
{
// Typedefs
// Hashtable to static class attributes addresses holding the
// runtime info switch values
typedef std::map<std::string, std::list<controlSwitches<int> *> >
ListInfoControlSwitches;

extern ListInfoControlSwitches* infoSwitchValues_;

class infoSwitch
    : public controlSwitches<int>
{
public:

    // Constructors

        //- Construct from name and value
        infoSwitch
        (
            const std::string& switchName,
            const int& switchValue,
            const std::string& switchDescription = ""
        )
        :
            controlSwitches<int>
            (
                switchName,
                debug::infoSwitchFromDict(switchName.c_str(), switchValue),
                switchDescription,
                debug::INFO_SWITCHES,
                &infoSwitchValues_
            )
        {}


    //- Destructor
    virtual ~infoSwitch()
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
