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

Description
    Write out all library debug switches

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "IFstream.H"
#include "IOobject.H"
#include "HashSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("new", "");
    argList::validOptions.insert("old", "");

    Foam::argList args(argc, argv);

    wordList currDebug(debug::debugSwitches().toc());
    wordList currInfo(debug::infoSwitches().toc());
    wordList currOpt(debug::optimisationSwitches().toc());
    wordList currTol(debug::tolerances().toc());
    wordList currConst(debug::constants().toc());

    if (args.optionFound("old") || args.optionFound("new"))
    {
        dictionary controlDict(IFstream(findEtcFile("controlDict", true))());

        wordHashSet oldDebug
        (
            controlDict.subDict("DebugSwitches").toc()
        );

        wordHashSet oldInfo
        (
            controlDict.subDict("InfoSwitches").toc()
        );

        wordHashSet oldOpt
        (
            controlDict.subDict("OptimisationSwitches").toc()
        );

        wordHashSet oldTol
        (
            controlDict.subDict("Tolerances").toc()
        );

        wordHashSet oldConst
        (
            controlDict.subDict("DimensionedConstants").toc()
        );


        wordHashSet hashset;
        wordList listing;


        // list old switches - but this can't work since the (old) inserted
        // switches are in both sets
        // Workaround:
        //  1. run without any options (get complete list)
        //  2. comment out DebugSwitches, run again with -new to find new ones
        //     and do a diff
        if (args.optionFound("old"))
        {
            IOobject::writeDivider(Info);

            hashset = wordHashSet(oldDebug);
            hashset -= wordHashSet(currDebug);
            listing = hashset.toc();
            sort(listing);
            Info<< "old DebugSwitches: " << listing << endl;

            hashset = wordHashSet(oldInfo);
            hashset -= wordHashSet(currInfo);
            listing = hashset.toc();
            sort(listing);
            Info<< "old InfoSwitches: " << listing << endl;

            hashset = wordHashSet(oldOpt);
            hashset -= wordHashSet(currOpt);
            listing = hashset.toc();
            sort(listing);
            Info<< "old OptimisationSwitches: " << listing << endl;

            hashset = wordHashSet(oldTol);
            hashset -= wordHashSet(currTol);
            listing = hashset.toc();
            sort(listing);
            Info<< "old Tolerances: " << listing << endl;

            hashset = wordHashSet(oldConst);
            hashset -= wordHashSet(currConst);
            listing = hashset.toc();
            sort(listing);
            Info<< "old DimensionedConstants: " << listing << endl;
        }

        // list new switches
        if (args.optionFound("new"))
        {
            IOobject::writeDivider(Info);

            hashset = wordHashSet(currDebug);
            hashset -= wordHashSet(oldDebug);

            listing = hashset.toc();
            sort(listing);
            Info<< "new DebugSwitches: " << listing << endl;

            hashset = wordHashSet(currInfo);
            hashset -= wordHashSet(oldInfo);
            listing = hashset.toc();
            sort(listing);
            Info<< "new InfoSwitches: " << listing << endl;

            hashset = wordHashSet(currOpt);
            hashset -= wordHashSet(oldOpt);
            listing = hashset.toc();
            sort(listing);
            Info<< "new OptimisationSwitches: " << listing << endl;

            hashset = wordHashSet(currTol);
            hashset -= wordHashSet(oldTol);
            listing = hashset.toc();
            sort(listing);
            Info<< "new Tolerances: " << listing << endl;

            hashset = wordHashSet(currConst);
            hashset -= wordHashSet(oldConst);
            listing = hashset.toc();
            sort(listing);
            Info<< "new DimensionedConstants: " << listing << endl;
        }
    }
    else
    {
        //IOobject::writeDivider(Info);

        sort(currDebug);

        Info << endl << "DebugSwitches: " << endl;
	forAll(currDebug, dI)
	{
	    Info << "  "
		 << currDebug[dI]
		 << " : "
		 << debug::debugSwitchFromDict(currDebug[dI].c_str(), 0)
		 << endl;
	}

        sort(currInfo);
        Info << endl << "InfoSwitches: " << endl;
	forAll(currInfo, iI)
	{
	    Info << "  "
		 << currInfo[iI]
		 << " : "
		 << debug::infoSwitchFromDict(currInfo[iI].c_str(), 0)
		 << endl;
	}

        sort(currOpt);
        Info << endl << "OptimisationSwitches: " << endl;
	forAll(currOpt, oI)
	{
	    if (currOpt[oI] == "commsType")
	    {
		token commsTypeValue;
		debug::optimisationSwitches().lookup("commsType", false, false).read(commsTypeValue);
	        Info << "  " << "commsType : " << commsTypeValue << endl;
	    }
	    else
	    {
		Info << "  "
		     << currOpt[oI]
		     << " : "
		     << debug::optimisationSwitchFromDict(currOpt[oI].c_str(), 0)
		     << endl;
	    }
	}

        sort(currTol);
        Info << endl << "Tolerances: " << endl;
	forAll(currTol, tI)
	{
	    Info << "  "
		 << currTol[tI]
		 << " : "
		 << debug::tolerancesFromDict(currTol[tI].c_str(), 0)
		 << endl;
	}

        sort(currConst);
        Info << endl << "Dimensioned Constants: " << endl;
	forAll(currConst, tI)
	{
	    Info << "  "
		 << currConst[tI]
		 << " : "
		 << debug::constantsFromDict(currConst[tI].c_str(), 0)
		 << endl;
	}
    }

    Info << endl << "Done." << endl;

    return 0;
}


// ************************************************************************* //
