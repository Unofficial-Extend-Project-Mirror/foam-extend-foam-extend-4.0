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

Namespace
    Foam::debug

Description
    Namespace for handling debugging switches.

SourceFiles
    debug.C

\*---------------------------------------------------------------------------*/

#ifndef debug_H
#define debug_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace debug
{
    //- The central control dictionary.
    //  Located in $WM_PROJECT_DIR/etc
    //  @sa Foam::findEtcFile()
    dictionary& controlDict();

    //- The DebugSwitches sub-dictionary in the central controlDict.
    dictionary& debugSwitches();

    //- The InfoSwitches sub-dictionary in the central controlDict.
    dictionary& infoSwitches();

    //- The OptimisationSwitches sub-dictionary in the central controlDict.
    dictionary& optimisationSwitches();

    //- The Tolerances sub-dictionary in the central controlDict.
    dictionary& tolerances();

    //- Lookup debug switch or add default value.
    int debugSwitch(const char* name, const int defaultValue = 0);

    //- Lookup info switch or add default value.
    int infoSwitch(const char* name, const int defaultValue = 0);

    //- Lookup optimisation switch or add default value.
    int optimisationSwitch(const char* name, const int defaultValue = 0);

    //- Lookup tolerances switch or add default value.
    double tolerances(const char* name, const double defaultValue = 0);

    //- Internal function to lookup a sub-dictionary from controlDict.
    dictionary& switchSet(const char* subDictName, dictionary*& subDictPtr);

} // End namespace debug


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
