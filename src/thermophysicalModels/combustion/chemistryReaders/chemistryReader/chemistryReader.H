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
#include "sutherlandTransport.H"
#include "specieThermo.H"
#include "janafThermo.H"
#include "perfectGas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class chemistryReader Declaration
\*---------------------------------------------------------------------------*/

class chemistryReader
{

public:

    // Public data types

        typedef sutherlandTransport<specieThermo<janafThermo<perfectGas> > >
            reactionThermo;

        typedef Reaction<reactionThermo> reaction;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        chemistryReader(const chemistryReader&);

        //- Disallow default bitwise assignment
        void operator=(const chemistryReader&);


public:

    //- Runtime type information
    TypeName("chemistryReader");


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


    // Destructor

        virtual ~chemistryReader()
        {}


    // Member Functions

        virtual const speciesTable& species() const = 0;
        virtual const HashPtrTable<reactionThermo>& sThermo() const = 0;
        virtual const SLPtrList<reaction>& reactions() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
