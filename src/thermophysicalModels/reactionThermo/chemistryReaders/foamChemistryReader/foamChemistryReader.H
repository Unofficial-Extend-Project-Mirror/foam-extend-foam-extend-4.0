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
    Foam::foamChemistryReader

Description
    Chemistry reader for FOAM format

SourceFiles
    foamChemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef foamChemistryReader_H
#define foamChemistryReader_H

#include "chemistryReader.H"
#include "fileName.H"
#include "typeInfo.H"
#include "HashPtrTable.H"
#include "SLPtrList.H"
#include "labelList.H"
#include "speciesTable.H"
#include "atomicWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class foamChemistry Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class foamChemistryReader
:
    public chemistryReader<ThermoType>
{
        //- Table of the thermodynamic data given in the foamChemistry file
        HashPtrTable<ThermoType> speciesThermo_;

        //- Table of species
        speciesTable speciesTable_;

        //- List of the reactions
        SLPtrList<Reaction<ThermoType> > reactions_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        foamChemistryReader(const foamChemistryReader&);

        //- Disallow default bitwise assignment
        void operator=(const foamChemistryReader&);


public:

    //- Runtime type information
    TypeName("foamChemistryReader");


    // Constructors

        //- Construct from foamChemistry and thermodynamics file names
        foamChemistryReader
        (
            const fileName& reactionsFileName,
            const fileName& thermoFileName
        );

        //- Construct by getting the foamChemistry and thermodynamics file names
        //  from dictionary
        foamChemistryReader(const dictionary& thermoDict);


    //- Destructor
    virtual ~foamChemistryReader()
    {}


    // Member functions

        //- Table of species
        const speciesTable& species() const
        {
            return speciesTable_;
        }

        //- Table of the thermodynamic data given in the foamChemistry file
        const HashPtrTable<ThermoType>& sThermo() const
        {
            return speciesThermo_;
        }

        //- List of the reactions
        const SLPtrList<Reaction<ThermoType> >& reactions() const
        {
            return reactions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "foamChemistryReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
