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

\*---------------------------------------------------------------------------*/

#include "argList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::argList::executable() const
{
    return executable_;
}


inline const Foam::fileName& Foam::argList::rootPath() const
{
    return rootPath_;
}


inline const Foam::fileName& Foam::argList::caseName() const
{
    return case_;
}


inline const Foam::fileName& Foam::argList::globalCaseName() const
{
    return globalCase_;
}


inline const Foam::ParRunControl& Foam::argList::parRunControl() const
{
    return parRunControl_;
}


inline Foam::fileName Foam::argList::path() const
{
    return rootPath()/caseName();
}


inline const Foam::stringList& Foam::argList::args() const
{
    return args_;
}


inline Foam::stringList& Foam::argList::args()
{
    return args_;
}


inline const Foam::string& Foam::argList::arg(const label index) const
{
    return args_[index];
}


inline Foam::label Foam::argList::size() const
{
    return args_.size();
}


inline const Foam::HashTable<Foam::string>& Foam::argList::options() const
{
    return options_;
}


inline Foam::HashTable<Foam::string>& Foam::argList::options()
{
    return options_;
}


inline const Foam::string& Foam::argList::option(const word& opt) const
{
    return options_[opt];
}


inline bool Foam::argList::optionFound(const word& opt) const
{
    return options_.found(opt);
}


inline Foam::IStringStream Foam::argList::optionLookup(const word& opt) const
{
    return IStringStream(options_[opt]);
}


// * * * * * * * * * * * * Template Specializations  * * * * * * * * * * * * //

namespace Foam
{
    // Template specialization for string
    template<>
    inline Foam::string
    Foam::argList::argRead<Foam::string>(const label index) const
    {
        return args_[index];
    }

    // Template specialization for word
    template<>
    inline Foam::word
    Foam::argList::argRead<Foam::word>(const label index) const
    {
        return args_[index];
    }

    // Template specialization for fileName
    template<>
    inline Foam::fileName
    Foam::argList::argRead<Foam::fileName>(const label index) const
    {
        return args_[index];
    }

    // Template specialization for string
    template<>
    inline Foam::string
    Foam::argList::optionRead<Foam::string>(const word& opt) const
    {
        return options_[opt];
    }

    // Template specialization for word
    template<>
    inline Foam::word
    Foam::argList::optionRead<Foam::word>(const word& opt) const
    {
        return options_[opt];
    }

    // Template specialization for fileName
    template<>
    inline Foam::fileName
    Foam::argList::optionRead<Foam::fileName>(const word& opt) const
    {
        return options_[opt];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline T Foam::argList::argRead(const label index) const
{
    T val;

    IStringStream(args_[index])() >> val;
    return val;
}


template<class T>
inline T Foam::argList::optionRead(const word& opt) const
{
    T val;

    optionLookup(opt)() >> val;
    return val;
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& opt,
    T& val
) const
{
    if (optionFound(opt))
    {
        val = optionRead<T>(opt);
        return true;
    }
    else
    {
        return false;
    }
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& opt,
    T& val,
    const T& deflt
) const
{
    if (optionReadIfPresent<T>(opt, val))
    {
        return true;
    }
    else
    {
        val = deflt;
        return false;
    }
}


template<class T>
inline T Foam::argList::optionLookupOrDefault
(
    const word& opt,
    const T& deflt
) const
{
    if (optionFound(opt))
    {
        return optionRead<T>(opt);
    }
    else
    {
        return deflt;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::string& Foam::argList::operator[](const label index) const
{
    return args_[index];
}


inline const Foam::string& Foam::argList::operator[](const word& opt) const
{
    return options_[opt];
}


// ************************************************************************* //
