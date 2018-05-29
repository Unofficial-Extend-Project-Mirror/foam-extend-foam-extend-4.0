/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    OutputControlDictionary

Description
    Host template class used to control input/output for a given type. Used to
    enable a combination of run-time selection using the TypeName macro (see
    typeInfo.H) and automatic read/write provided by regIOobject part of the
    IOdictionary.

    For example of usage, see $FOAM_SRC/ODE/sixDOFODE/sixDOFODE.H

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "OutputControlDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class PolicyType>
Foam::OutputControlDictionary<PolicyType>::OutputControlDictionary
(
    const IOobject& io,
    const PolicyType& pt
)
:
    IOdictionary(io),
    pt_(pt)
{
    // Note: parameter pt may be incomplete here, must not call its member
    // functions inside the constructor body.
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class PolicyType>
Foam::OutputControlDictionary<PolicyType>::~OutputControlDictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class PolicyType>
bool Foam::OutputControlDictionary<PolicyType>::writeData(Ostream& os) const
{
    return pt_.writeData(os);
}


// ************************************************************************* //
