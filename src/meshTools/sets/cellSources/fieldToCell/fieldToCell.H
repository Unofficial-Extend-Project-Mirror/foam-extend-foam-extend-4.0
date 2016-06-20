/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::fieldToCell

Description
    A topoSetSource to select cells based on field values.

SourceFiles
    fieldToCell.C

\*---------------------------------------------------------------------------*/

#ifndef fieldToCell_H
#define fieldToCell_H

#include "topoSetSource.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fieldToCell Declaration
\*---------------------------------------------------------------------------*/

class fieldToCell
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;


        //- Name of volScalarField
        word fieldName_;

        //- min to subset field values with
        scalar min_;

        //- max ,,
        scalar max_;


    // Private Member Functions

        //- Depending on field values add to or delete from cellSet.
        void applyToSet
        (
            const topoSetSource::setAction action,
            const scalarField& field,
            topoSet& set
        ) const;

public:

    //- Runtime type information
    TypeName("fieldToCell");

    // Constructors

        //- Construct from components
        fieldToCell
        (
            const polyMesh& mesh,
            const word& fieldName,
            const scalar min,
            const scalar max
        );

        //- Construct from dictionary
        fieldToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        fieldToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    // Destructor

        virtual ~fieldToCell();


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
