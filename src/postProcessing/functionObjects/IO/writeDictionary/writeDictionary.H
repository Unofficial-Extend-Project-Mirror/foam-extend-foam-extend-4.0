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
    Foam::writeDictionary

Group
    grpUtilitiesFunctionObjects

Description
    This function object writes dictionaries on start-up, and on change

SourceFiles
    writeDictionary.C
    IOwriteDictionary.H

\*---------------------------------------------------------------------------*/

#ifndef writeDictionary_H
#define writeDictionary_H

#include "wordList.H"
#include "runTimeSelectionTables.H"
#include "SHA1Digest.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class writeDictionary Declaration
\*---------------------------------------------------------------------------*/

class writeDictionary
{
protected:

    // Private data

        //- Name of this set of writeDictionary
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- Names of dictionaries to monitor
        wordList dictNames_;

        //- List of changed dictionaries (only those registered to database)
        List<SHA1Digest> digests_;


    // Private Member Functions

        //- Helper function to write the dictionary if found at location
        bool tryDirectory
        (
            const label dictI,
            const word& location,
            bool& firstDict
        );

        //- Disallow default bitwise copy construct
        writeDictionary(const writeDictionary&);

        //- Disallow default bitwise assignment
        void operator=(const writeDictionary&);


public:

    //- Runtime type information
    TypeName("writeDictionary");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        writeDictionary
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~writeDictionary();


    // Member Functions

        //- Return name of the writeDictionary
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the writeDictionary data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the writeDictionary
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
