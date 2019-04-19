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
    Foam::dlLibraryTable

Description
    A table of dynamically loaded libraries

SourceFiles
    dlLibraryTable.C

\*---------------------------------------------------------------------------*/

#ifndef dlLibraryTable_H
#define dlLibraryTable_H

#include "label.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dlLibraryTable Declaration
\*---------------------------------------------------------------------------*/

class dlLibraryTable
{
    // Private Member Functions

        DynamicList<void*> libPtrs_;

        DynamicList<fileName> libNames_;


        //- Disallow default bitwise copy construct
        dlLibraryTable(const dlLibraryTable&);

        //- Disallow default bitwise assignment
        void operator=(const dlLibraryTable&);


public:

    // Declare name of the class and its debug switch
    ClassName("dlLibraryTable");

    // Constructors

        //- Construct null
        dlLibraryTable();

        //- Construct from dictionary and name of 'libs' entry giving
        //  the libraries to load
        dlLibraryTable(const dictionary&, const word&);


    //- Destructor
    ~dlLibraryTable();


    // Member Functions

        //- Open the named library, optionally with warnings if problems occur
        bool open(const fileName& name, const bool verbose = true);

        //- Close the named library, optionally with warnings if problems occur
        bool close(const fileName& name, const bool verbose = true);

        //- Find the handle of the named library
        void* findLibrary(const fileName& name);

        //- Open all the libraries listed in the 'libsEntry' entry in the
        //  given dictionary if present
        bool open(const dictionary&, const word& libsEntry);

        //- Open all the libraries listed in the 'libsEntry' entry in the
        //  given dictionary if present and check the additions
        //  to the given constructor table
        template<class TablePtr>
        bool open
        (
            const dictionary&,
            const word& libsEntry,
            const TablePtr& tablePtr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dlLibraryTableTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
