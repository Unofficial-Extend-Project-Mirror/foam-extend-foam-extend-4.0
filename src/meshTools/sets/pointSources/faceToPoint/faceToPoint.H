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
    Foam::faceToPoint

Description
    A topoSetSource to select points based on usage in faces.

SourceFiles
    faceToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef faceToPoint_H
#define faceToPoint_H

#include "topoSetSource.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceToPoint Declaration
\*---------------------------------------------------------------------------*/

class faceToPoint
:
    public topoSetSource
{

public:
        //- Enumeration defining the valid options
        enum faceAction
        {
            ALL
        };

private:

        //- Add usage string
        static addToUsageTable usage_;

        static const NamedEnum<faceAction, 1> faceActionNames_;

        //- Name of set to use
        word setName_;

        //- Option
        faceAction option_;


    // Private Member Functions

        //- Depending on face to cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("faceToPoint");

    // Constructors

        //- Construct from components
        faceToPoint
        (
            const polyMesh& mesh,
            const word& setName,
            const faceAction option
        );

        //- Construct from dictionary
        faceToPoint
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        faceToPoint
        (
            const polyMesh& mesh,
            Istream&
        );


    // Destructor

        virtual ~faceToPoint();


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
