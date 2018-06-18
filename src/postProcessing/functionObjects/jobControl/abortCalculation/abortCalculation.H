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
    Foam::abortCalculation

Group
    grpJobControlFunctionObjects

Description
    Watches for presence of the named file in the $FOAM_CASE directory
    and aborts the calculation if it is present.

    Currently the following action types are supported:
    - noWriteNow
    - writeNow
    - nextWrite

SourceFiles
    abortCalculation.C
    IOabortCalculation.H

\*---------------------------------------------------------------------------*/

#ifndef abortCalculation_H
#define abortCalculation_H

#include "NamedEnum.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                      Class abortCalculation Declaration
\*---------------------------------------------------------------------------*/

class abortCalculation
{
public:

    // Public data

        //- Enumeration defining the type of action
        enum actionType
        {
            noWriteNow,    /*!< stop immediately without writing data */
            writeNow,      /*!< write data and stop immediately */
            nextWrite      /*!< stop the next time data are written */
        };

private:

    // Private data

        //- Name of the abort file unless otherwise specified
        word name_;

        const objectRegistry& obr_;

        //- The fully-qualified name of the abort file
        fileName abortFile_;

        //- Action type names
        static const NamedEnum<actionType, 3> actionTypeNames_;

        //- The type of action
        actionType action_;


    // Private Member Functions

        //- Remove abort file.
        void removeFile() const;

        //- Disallow default bitwise copy construct
        abortCalculation(const abortCalculation&);

        //- Disallow default bitwise assignment
        void operator=(const abortCalculation&);


public:

    //- Runtime type information
    TypeName("abort");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        abortCalculation
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFilesUnused = false
        );


    //- Destructor
    virtual ~abortCalculation();


    // Member Functions

        //- Return name of the abort file
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the dictionary settings
        virtual void read(const dictionary&);

        //- Execute, check existence of abort file and take action
        virtual void execute();

        //- Execute at the final time-loop, used for cleanup
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Execute, check existence of abort file and take action
        virtual void write();

        //- Update for changes of mesh - does nothing
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh - does nothing
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
