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

\*---------------------------------------------------------------------------*/

#include "regExp.H"
#include "label.H"
#include "foamString.H"
#include "List.H"
#include "IOstreams.H"

// Alternative regular expression libraries to consider are:
// Boost http://www.boost.org/libs/regex/doc/
// GRETA http://research.microsoft.com/projects/greta/
// Henry Spencer's http://arglist.com/regex/
//
// Chose DEELX http://www.regexlab.com/en/deelx/
// for its ease of integration - one header file
#include "deelx.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regExp::regExp()
:
    preg_(0)
{}


Foam::regExp::regExp(const char* pattern, const bool ignoreCase)
:
    preg_(0)
{
    set(pattern, ignoreCase);
}


Foam::regExp::regExp(const std::string& pattern, const bool ignoreCase)
:
    preg_(0)
{
    set(pattern.c_str(), ignoreCase);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regExp::~regExp()
{
    clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regExp::set(const char* pattern, const bool ignoreCase) const
{
    clear();

    // avoid NULL pointer and zero-length patterns
    if (pattern && *pattern)
    {
        int cflags = EXTENDED;
        if (ignoreCase)
        {
            cflags |= IGNORECASE;
        }

        preg_ = new CRegexpT<char>(pattern, cflags);
    }
}


void Foam::regExp::set(const std::string& pattern, const bool ignoreCase) const
{
    return set(pattern.c_str(), ignoreCase);
}


bool Foam::regExp::clear() const
{
    if (preg_)
    {
        delete preg_;
        preg_ = 0;

        return true;
    }

    return false;
}


std::string::size_type Foam::regExp::find(const std::string& str) const
{
    std::string::size_type pos = std::string::npos;

    if (preg_ && !str.empty())
    {
        const MatchResult result = preg_->Match(str.c_str());

        if( 0 < result.IsMatched() )
        {
            pos = result.GetStart();
        }
    }

    return pos;
}


bool Foam::regExp::match(const std::string& str) const
{
    bool isExactMatch = false;


    if( preg_ && !str.empty() )
    {
        const MatchResult result = preg_->MatchExact(str.c_str());
        isExactMatch = (0 < result.IsMatched());
    }

    return isExactMatch;
}


bool Foam::regExp::match(const string& str, List<string>& groups) const
{
    bool isMatch = false;

    if( preg_ && !str.empty() )
    {
        const MatchResult results = preg_->MatchExact(str.c_str());
        isMatch = (0 < results.IsMatched());

        if( isMatch )
        {
            int const notFound = -1;
            int start, end;
            const int groupsCount = results.MaxGroupNumber();
            groups.setSize(groupsCount);

            for (int i = 0; groupsCount > i; ++i)
            {
                start = results.GetGroupStart(i);
                end   = results.GetGroupEnd(i);

                if ((notFound < start) && (notFound < end))
                {
                    groups[i] = str.substr(start, end - start);
                }
                else
                {
                    groups[i].clear();
                }
            }
        }
    }

    if( !isMatch )
    {
        groups.clear();
    }

    return isMatch;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::regExp::operator=(const char* pat)
{
    set(pat);
}


void Foam::regExp::operator=(const std::string& pat)
{
    set(pat);
}


// ************************************************************************* //
