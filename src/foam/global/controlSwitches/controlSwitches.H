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
    Foam::controlSwitches

Author
    Martin Beaudoin, Hydro-Quebec, 2014.  All rights reserved

Description
    Handling of runtime controlSwitches
    - Debug switches
    - Info switches
    - Optimisation switches
    - Tolerances switches
    - DimensionedConstants values

    This class is a specialization of the safeBool class, where additional
    boolean operators are allowed. This class is also offering access to
    information like the switch name, etc.

Note:
    We will be using std::map and std::list for the containers for this class.
    In some sections of the code where controlSwitches are being used, we are
    at a relatively low level, so we cannot use HashTable and DynamicList.
    MB 23/05/2014

SourceFiles
    controlSwitchesI.H
    controlSwitches.C
    controlSwitchesIO.C

\*---------------------------------------------------------------------------*/

#ifndef controlSwitches_H
#define controlSwitches_H

#include <list>
#include <map>
#include <string>
#include "safeBool.H"

#include "scalar.H"
#include "debug.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace debug
{
    template<class T> class controlSwitches;
}

template<class T>
inline bool operator==
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator==
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator==
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
);

template<class T>
inline bool operator!=
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator!=
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator!=
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
);

template<class T>
inline bool operator>=
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator>=
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator>=
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
);

template<class T>
inline bool operator>
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator>
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator>
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
);

template<class T>
inline bool operator<=
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator<=
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator<=
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
);

template<class T>
inline bool operator<
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator<
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
);

template<class T>
inline bool operator<
(
    const debug::controlSwitches<T>&lhs,
    const T& rhs
);

namespace debug
{

/*---------------------------------------------------------------------------*\
                        Class controlSwitches Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class controlSwitches
:
    public safeBool<controlSwitches<T> > // using the CRTP idiom
{
    // Private data

        //- Switch section
        globalControlDictSwitchSet switchSet_;

        //- Switch name
        std::string switchName_;

        //- Switch value
        T switchValue_;

        //- Switch Description
        std::string switchDescription_;

        //- Handle to runTime switches list
        std::map<std::string, std::list<controlSwitches<T> *> >*
            switchValuesTable_;


    // Private member functions

        //- Disallow construct as copy
        controlSwitches(const controlSwitches<T>&);

        //- Disallow default bitwise assignment
        void operator=(const controlSwitches<T>&);


public:

    // Constructors

        controlSwitches
        (
            const std::string& switchName,
            const T& switchValue,
            const std::string& switchDescription,
            globalControlDictSwitchSet switchSet,
            std::map<std::string, std::list<controlSwitches<T> *> >**
            switchesValues
        );


    //- Destructor
    ~controlSwitches();


    // Member Functions

        // Check
        bool test() const
        {
            // Perform Boolean logic here
            return switchValue_ != 0;
        }


    // Member Operators

        //- Assignement operator
        void operator=(const T&);

        //- & operator
        const T operator&(const T&);

        //- Return value
        inline T operator()() const
        {
            return switchValue_;
        }

        // Return switch name
        inline std::string switchName()
        {
            return switchName_;
        }

        //- Return switch description
        inline std::string switchDescription()
        {
            return switchDescription_;
        }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Print the content of a controlSwitches list to the console.
// We are forced to pass the string descriptions of the Pstream::commsTypes
// when dealing with the optimisationSwitches group because this switch is in
// fact an enum but we need to specify its corresponding string equivalent in
// a controlDict dictionary. And at the low level we are playing, including
// Pstream.H is out of the question.
// MB 2015
template<class T>
void printControlSwitches
(
    const std::string&,
    const std::map<std::string, std::list<debug::controlSwitches<T> *> >*,
    const char* commsTypesNames[] = NULL

);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace debug


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "controlSwitchesI.H"

#ifdef NoRepository
#   include "controlSwitches.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
