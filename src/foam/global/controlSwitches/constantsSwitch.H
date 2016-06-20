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
    Specialization of controlSwitches class for dimensionedConstant switches

\*---------------------------------------------------------------------------*/

#ifndef constantsSwitch_H
#define constantsSwitch_H

#include "controlSwitches.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace debug
{

// Typedefs
// Hashtables to static class attributes addresses holding the
// runtime dimensionedConstants values
typedef std::map<std::string, std::list<controlSwitches<Foam::scalar> *> >
ListConstantsControlSwitches;

extern ListConstantsControlSwitches* constantsSwitchValues_;


class constantsSwitch
:
    public controlSwitches<scalar>
{
    // Private member functions

        //- Disallow construct as copy
        constantsSwitch(const constantsSwitch&);

        //- Disallow default bitwise assignment
        void operator=(const constantsSwitch&);


public:

    // Constructors

        //- Construct from name and value
        constantsSwitch
        (
            const std::string& switchName,
            const Foam::scalar& switchValue,
            const std::string& switchDescription = ""
        )
        :
            controlSwitches<scalar>
            (
                switchName,
                debug::constantsFromDict(switchName.c_str(), switchValue),
                switchDescription,
                debug::DIMENSIONED_CONSTANTS,
                &constantsSwitchValues_
            )
        {}


    //- Destructor
    virtual ~constantsSwitch()
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
