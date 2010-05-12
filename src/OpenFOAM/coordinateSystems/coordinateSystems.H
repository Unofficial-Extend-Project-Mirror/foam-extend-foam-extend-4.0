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

Class
    Foam::coordinateSystems

Description
    Provides a centralized coordinateSystem collection.

SourceFiles
    coordinateSystems.C

\*---------------------------------------------------------------------------*/
#ifndef coordinateSystems_H
#define coordinateSystems_H

#include "coordinateSystem.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateSystems Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystems
:
    public PtrList<coordinateSystem>
{

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const coordinateSystems&);

public:

    // Static data members

        //- Type name of list contents
        static const word dataType;

    // Public Member Functions


    // Constructors

        //- Construct null
        coordinateSystems();

        //- Read construct from registry, name. instance
        coordinateSystems
        (
            const objectRegistry& registry,
            const word& name = "coordinateSystems",
            const fileName& instance = "constant"
        );

        //- Read construct from IOobject
        coordinateSystems(const IOobject&);

    // Member Functions


        //- Find and return index for a given keyword,  returns -1 if not found
        label find(const word& key) const;

        //- Search for given keyword
        bool found(const word& keyword) const;

        //- Return the table of contents (list of all keywords)
        wordList toc() const;

        //- Rewrite coordinateSystem entry with appropriate dictionary entry
        //
        //  This replaces coordinateSystems entries with their contents.
        //  For example,
        //  @verbatim
        //  cat1
        //  {
        //      coordinateSystem    system_10;
        //      porosity        0.781;
        //      Darcy
        //      {
        //          d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
        //          f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
        //      }
        //  }
        //  @endverbatim
        //  might get re-written as the following (depending on the value of
        //  @c system_10 in the @c constant/coordinateSystems file):
        //  @verbatim
        //  cat1
        //  {
        //      coordinateSystem
        //      {
        //          origin  (0 0 0);
        //          e3      (1 0 0);
        //          e1      (0 0 -1);
        //      }
        //      porosity        0.781;
        //      Darcy
        //      {
        //          d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
        //          f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
        //      }
        //  }
        //  @endverbatim
        //  When this form of re-writing is used, the coordinateRotation is
        //  reduced to the axes specification.
        bool rewriteDict(dictionary& dict, bool noType = false) const;

        //- write data
        bool writeData(Ostream&, bool subDict = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
