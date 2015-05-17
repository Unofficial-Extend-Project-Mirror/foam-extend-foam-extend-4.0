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

Class
    Foam::chemistryReader

Description
    Abstract class for reading chemistry

SourceFiles
    chemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryReader_H
#define chemistryReader_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class chemistryReader Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryReader
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        chemistryReader(const chemistryReader&);

        //- Disallow default bitwise assignment
        void operator=(const chemistryReader&);


public:

    //- Runtime type information
    TypeName("chemistryReader");

    //- The type of thermo package the reader was instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct null
        chemistryReader()
        {}


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistryReader,
            dictionary,
            (
                const dictionary& thermoDict
            ),
            (thermoDict)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<chemistryReader> New(const dictionary& thermoDict);


    //- Destructor
    virtual ~chemistryReader()
    {}


    // Member Functions

        virtual const speciesTable& species() const = 0;

        // Bug fix: cannot be named speciesThermo because of name clash
        // with the speciesThermo class.  HJ, 21/Sep/2010
        virtual const HashPtrTable<ThermoType>& sThermo() const = 0;

        virtual const SLPtrList<Reaction<ThermoType> >& reactions() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistryReader(Thermo)                                           \
                                                                              \
defineTemplateTypeNameAndDebug                                                \
(                                                                             \
    chemistryReader<Thermo>,                                                  \
    0                                                                         \
);                                                                            \
                                                                              \
defineTemplateRunTimeSelectionTable(chemistryReader<Thermo>, dictionary);


#define makeChemistryReaderType(SS, Thermo)                                   \
                                                                              \
defineNamedTemplateTypeNameAndDebug(SS<Thermo>, 0);                           \
                                                                              \
chemistryReader<Thermo>::adddictionaryConstructorToTable<SS<Thermo> >         \
    add##SS##Thermo##ConstructorToTable_;


#define addChemistryReaderType(SS, Thermo)                                    \
                                                                              \
defineTypeNameAndDebug(SS, 0);                                                \
                                                                              \
chemistryReader<Thermo>::adddictionaryConstructorToTable<SS>                  \
    add##SS##Thermo##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistryReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
