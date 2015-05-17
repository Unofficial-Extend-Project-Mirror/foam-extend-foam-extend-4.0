/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "dummyControlDict.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dummyControlDict, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dummyControlDict::dummyControlDict()
{
    this->set("deltaT", 1);
    this->set("writeFrequency", 1);
}


Foam::dummyControlDict::dummyControlDict(const fileName& mcdFile)
{
    this->set("deltaT", 1);
    this->set("writeFrequency", 1);

    IFstream mcdStream(mcdFile);
    if (!mcdStream.good())
    {
        FatalErrorIn("dummyControlDict::dummyControlDict")
            << "Cannot find the multiControlDict file " << mcdFile << "."
            << abort(FatalError);
    }
    dictionary mcdDict(mcdStream);
    if (mcdDict.subDict("multiSolverControl").found("timeFormat"))
    {
        word tf(mcdDict.subDict("multiSolverControl").lookup("timeFormat"));
        this->set("timeFormat", tf);
    }
    if (mcdDict.subDict("multiSolverControl").found("timePrecision"))
    {
        label tp
        (
            readLabel
            (
                mcdDict.subDict("multiSolverControl").lookup("timePrecision")
            )
        );
        this->set("timePrecision", tp);
    }
}


Foam::dummyControlDict::dummyControlDict(const dictionary& mcdDict)
{
    this->set("deltaT", 1);
    this->set("writeFrequency", 1);

    if (mcdDict.subDict("multiSolverControl").found("timeFormat"))
    {
        word tf(mcdDict.subDict("multiSolverControl").lookup("timeFormat"));
        this->set("timeFormat", tf);
    }
    if (mcdDict.subDict("multiSolverControl").found("timePrecision"))
    {
        label tp
        (
            readLabel
            (
                mcdDict.subDict("multiSolverControl").lookup("timePrecision")
            )
        );
        this->set("timePrecision", tp);
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dummyControlDict::~dummyControlDict()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
