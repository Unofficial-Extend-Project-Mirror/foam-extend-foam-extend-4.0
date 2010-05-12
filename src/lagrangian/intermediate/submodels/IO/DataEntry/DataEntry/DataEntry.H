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
    Foam::DataEntry

Description

SourceFiles
    DataEntry.C
    NewDataEntry.C

\*---------------------------------------------------------------------------*/

#ifndef DataEntry_H
#define DataEntry_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class DataEntry Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DataEntry
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        DataEntry(const DataEntry<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const DataEntry<Type>&);


protected:

    // Protected data

        //- Coefficients dictionary
        const dictionary dict_;

        //- Name of entry
        const word entry_;


public:

    //- Runtime type information
    TypeName("DataEntry")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DataEntry,
        dictionary,
        (
            const word& entryName,
            const dictionary& dict
        ),
        (entryName, dict)
    );


    // Constructor

        //- Construct from type name and dictionary
        DataEntry
        (
            const word& TypeName,
            const word& entryName,
            const dictionary& dict
        );


    //- Selector
        static autoPtr<DataEntry<Type> > New
        (
            const word& entryName,
            const dictionary& dict
        );


    //- Destructor

        virtual ~DataEntry();


    // Member Functions

        // Access

            //- Return the dictionary
            const dictionary& dict() const;


        // Evaluation

            //- Return value as a function of (scalar) independent variable
            virtual Type value(const scalar x) const = 0;

            //- Integrate between two (scalar) values
            virtual Type integrate(const scalar x1, const scalar x2) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDataEntry(Type)                                                   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(DataEntry<Type>, 0);                  \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        DataEntry<Type>,                                                      \
        dictionary                                                            \
    );


#define makeDataEntryType(SS, Type)                                           \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                         \
                                                                              \
    DataEntry<Type>::adddictionaryConstructorToTable<SS<Type> >               \
        add##SS##Type##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DataEntry.C"
#   include "NewDataEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
