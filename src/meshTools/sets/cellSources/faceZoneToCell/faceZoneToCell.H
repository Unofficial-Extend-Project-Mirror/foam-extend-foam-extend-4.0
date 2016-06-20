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
    Foam::faceZoneToCell

Description
    A topoSetSource to select cells based on side of faceZone.

SourceFiles
    faceZoneToCell.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneToCell_H
#define faceZoneToCell_H

#include "topoSetSource.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceZoneToCell Declaration
\*---------------------------------------------------------------------------*/

class faceZoneToCell
:
    public topoSetSource
{
public:
        //- Enumeration defining the valid options
        enum faceAction
        {
            MASTER,
            SLAVE
        };

private:

    // Private data

        static const NamedEnum<faceAction, 2> faceActionNames_;

        //- Add usage string
        static addToUsageTable usage_;

        //- Name/regular expression of faceZone
        wordRe zoneName_;

        //- Option
        faceAction option_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("faceZoneToCell");

    // Constructors

        //- Construct from components
        faceZoneToCell
        (
            const polyMesh& mesh,
            const word& zoneName,
            const faceAction option
        );

        //- Construct from dictionary
        faceZoneToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        faceZoneToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    // Destructor

        virtual ~faceZoneToCell();


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
