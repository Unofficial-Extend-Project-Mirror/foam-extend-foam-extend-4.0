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
    Foam::referredCellList

Description

SourceFiles
    referredCellListI.H
    referredCellList.C

\*---------------------------------------------------------------------------*/

#ifndef referredCellList_H
#define referredCellList_H

#include "referredCell.H"
#include "molecule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class interactionLists;

/*---------------------------------------------------------------------------*\
                      Class referredCellList Declaration
\*---------------------------------------------------------------------------*/

class referredCellList
:
    public List<referredCell>
{
    // Private data

        const interactionLists& il_;


    // Private Member Functions

        void buildReferredCellList
        (
            bool pointPointListBuild
        );


public:

    // Constructors

        //- Construct lists by searching the mesh
        referredCellList
        (
            interactionLists& il,
            bool pointPointListBuild
        );

        //- Construct from file
        referredCellList (interactionLists& il);


    // Destructor

        ~referredCellList();


    // Member Functions

        void referMolecules(const List<DynamicList<molecule*> >& cellOccupancy);

        inline const interactionLists& il() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "referredCellListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
