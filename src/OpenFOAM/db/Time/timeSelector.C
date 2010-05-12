/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "timeSelector.H"
#include "ListOps.H"
#include "argList.H"
#include "Time.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeSelector::timeSelector()
:
    scalarRanges()
{}


Foam::timeSelector::timeSelector(Istream& is)
:
    scalarRanges(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::timeSelector::selected(const instant& value) const
{
    return scalarRanges::selected(value.value());
}


Foam::List<bool> Foam::timeSelector::selected(const List<instant>& Times) const
{
    List<bool> lst(Times.size(), false);

    // check ranges
    forAll(Times, i)
    {
        if (selected(Times[i]))
        {
            lst[i] = true;
        }
    }

    // avoid false positive on "constant"
    forAll(Times, i)
    {
        if (Times[i].name() == "constant")
        {
            lst[i] = false;
            break;
        }
    }

    // check specific values
    forAll(*this, rangeI)
    {
        if (operator[](rangeI).isExact())
        {
            scalar target = operator[](rangeI).value();

            int nearestIndex = -1;
            scalar nearestDiff = Foam::GREAT;

            forAll(Times, timeIndex)
            {
                if (Times[timeIndex].name() == "constant") continue;

                scalar diff = fabs(Times[timeIndex].value() - target);
                if (diff < nearestDiff)
                {
                    nearestDiff = diff;
                    nearestIndex = timeIndex;
                }
            }

            if (nearestIndex >= 0)
            {
                lst[nearestIndex] = true;
            }
        }
    }

    return lst;
}


Foam::List<Foam::instant> Foam::timeSelector::select
(
    const List<instant>& Times
) const
{
    return subset(selected(Times), true, Times);
}


void Foam::timeSelector::inplaceSelect
(
    List<instant>& Times
) const
{
    inplaceSubset(selected(Times), true, Times);
}


void Foam::timeSelector::addOptions
(
    const bool constant,
    const bool zeroTime
)
{
    if (constant)
    {
        argList::validOptions.insert("constant", "");
    }
    if (zeroTime)
    {
        argList::validOptions.insert("zeroTime", "");
    }
    argList::validOptions.insert("noZero", "");
    argList::validOptions.insert("time", "time");
    argList::validOptions.insert("latestTime", "");
}


Foam::List<Foam::instant> Foam::timeSelector::select
(
    const List<instant>& timeDirs,
    const argList& args
)
{
    if (timeDirs.size())
    {
        List<bool> selectTimes(timeDirs.size(), true);

        if (args.options().found("time"))
        {
            selectTimes = timeSelector
            (
                IStringStream(args.options()["time"])()
            ).selected(timeDirs);
        }
        else if (args.options().found("latestTime"))
        {
            selectTimes = false;

            // avoid false match on constant/ or 0/
            if (timeDirs.size() > 2)
            {
                selectTimes[timeDirs.size() - 1] = true;
            }
        }

        if (timeDirs.size() > 1)
        {
            selectTimes[0] = args.options().found("constant");

            if (args.options().found("noZero"))
            {
                selectTimes[1] = false;
            }
            else if (argList::validOptions.found("zeroTime"))
            {
                selectTimes[1] = args.options().found("zeroTime");
            }
        }

        return subset(selectTimes, true, timeDirs);
    }
    else
    {
        return timeDirs;
    }
}


Foam::List<Foam::instant> Foam::timeSelector::select0
(
    Time& runTime,
    const argList& args
)
{
    instantList timeDirs = timeSelector::select(runTime.times(), args);

    if (!timeDirs.size())
    {
        FatalErrorIn(args.executable())
            << "No times selected"
            << exit(FatalError);
    }

    runTime.setTime(timeDirs[0], 0);

    return timeDirs;
}


// ************************************************************************* //
