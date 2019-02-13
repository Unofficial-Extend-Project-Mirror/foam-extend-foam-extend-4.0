/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    sortEdgesIntoChains

Description
    Creates boundary faces from a chain of open edges

SourceFiles
    sortEdgesIntoChains.C

\*---------------------------------------------------------------------------*/

#ifndef sortEdgesIntoChains_H
#define sortEdgesIntoChains_H

#include "DynList.H"
#include "labelList.H"
#include "edge.H"
#include "Map.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sortEdgesIntoChains Declaration
\*---------------------------------------------------------------------------*/

class sortEdgesIntoChains
{
    // Members
        const DynList<edge>& bEdges_;

        bool openEdges_;

        Map<label> newNodeLabel_;

        DynList<DynList<label> > edgesAtPoint_;

        DynList<DynList<label> > createdChains_;

    // Private member functions
        void createNodeLabels();

        bool findPointsBelongingToTheChain
        (
            const label currPos,
            DynList<bool>& chainEdges
        ) const;

        void shrinkEdges(const DynList<bool>& chainEdges);

        void createChainFromEdges(const DynList<bool>& chainEdges);

        void sortEdges();

    public:

        sortEdgesIntoChains
        (
            const DynList<edge>& bEdges
        );

        ~sortEdgesIntoChains();

    // Member functions

        //- a list of points which have not yet been resolved
        const DynList<DynList<label> >& sortedChains() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
