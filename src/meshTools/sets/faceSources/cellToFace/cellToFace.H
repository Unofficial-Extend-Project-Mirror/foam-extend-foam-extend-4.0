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
    Foam::cellToFace

Description
    A topoSetSource to select a faceSet from cells.

    Either all faces of cell or some other criterion.
    See implementation.
    Note: when picking up coupled faces uses cells on neighbouring processors.

SourceFiles
    cellToFace.C

\*---------------------------------------------------------------------------*/

#ifndef cellToFace_H
#define cellToFace_H

#include "topoSetSource.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellToFace Declaration
\*---------------------------------------------------------------------------*/

class cellToFace
:
    public topoSetSource
{
public:
        //- Enumeration defining the valid options
        enum cellAction
        {
            ALL,
            BOTH
        };


private:

        //- Add usage string
        static addToUsageTable usage_;

        static const NamedEnum<cellAction, 2> cellActionNames_;

        //- Name of set to use
        word setName_;

        //- Option
        cellAction option_;


    // Private Member Functions

        //- Depending on face to cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("cellToFace");

    // Constructors

        //- Construct from components
        cellToFace
        (
            const polyMesh& mesh,
            const word& setName,
            const cellAction option
        );

        //- Construct from dictionary
        cellToFace
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        cellToFace
        (
            const polyMesh& mesh,
            Istream&
        );


    // Destructor

        virtual ~cellToFace();


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
