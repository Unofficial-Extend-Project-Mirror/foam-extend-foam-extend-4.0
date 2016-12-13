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
    Foam::infoSwitch

Author
    Martin Beaudoin, Hydro-Quebec, 2014.  All rights reserved

Description
    Handling of runtime controlSwitches.
    Specialization of controlSwitches class for tolerances switches

\*---------------------------------------------------------------------------*/

#ifndef tolerancesSwitch_H
#define tolerancesSwitch_H

#include "controlSwitches.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace debug
{

// Typedefs
// Hashtable to static class attributes addresses holding the
// runtime tolerances switch values
typedef std::map<std::string, std::list<controlSwitches<Foam::scalar> *> >
ListTolerancesControlSwitches;

extern ListTolerancesControlSwitches* tolerancesSwitchValues_;

class tolerancesSwitch
:
    public controlSwitches<scalar>
{
    // Private member functions

        //- Disallow construct as copy
        tolerancesSwitch(const tolerancesSwitch&);

        //- Disallow default bitwise assignment
        void operator=(const tolerancesSwitch&);


public:

    // Constructors

        //- Construct from name and value
        tolerancesSwitch
        (
            const std::string& switchName,
            const Foam::scalar& switchValue,
            const std::string& switchDescription = ""
        )
        :
            controlSwitches<scalar>
            (
                switchName,
                debug::tolerancesFromDict(switchName.c_str(), switchValue),
                switchDescription,
                debug::TOLERANCES,
                &tolerancesSwitchValues_
            )
        {}


    //- Destructor
    virtual ~tolerancesSwitch()
    {}


    // Member operators

        //- Assignement operator
        void operator=(const Foam::scalar& rhs)
        {
            controlSwitches<Foam::scalar>::operator=(rhs);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace debug


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
