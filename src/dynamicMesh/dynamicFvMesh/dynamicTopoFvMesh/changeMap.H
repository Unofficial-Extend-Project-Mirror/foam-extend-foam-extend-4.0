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
    changeMap

Description
    Accumulate topology change statistics

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    changeMapI.H

\*---------------------------------------------------------------------------*/

#ifndef changeMap_H
#define changeMap_H

#include "objectMap.H"
#include "labelList.H"
#include "dictionary.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class changeMap;

// * * * * * * * * Forward declaration of friend fuctions  * * * * * * * * * //

Ostream& operator<<(Ostream&, const changeMap&);


/*---------------------------------------------------------------------------*\
                           Class changeMap Declaration
\*---------------------------------------------------------------------------*/

class changeMap
:
    public dictionary
{
    // Entity index
    label index_;

    // Coupled patch index
    label pIndex_;

    // Type
    label type_;

    // Entities that were added during the operation.
    DynamicList<objectMap> addedPoints_;
    DynamicList<objectMap> addedEdges_;
    DynamicList<objectMap> addedFaces_;
    DynamicList<objectMap> addedCells_;

    // Entities that were removed during the operation
    DynamicList<label> removedPoints_;
    DynamicList<label> removedEdges_;
    DynamicList<label> removedFaces_;
    DynamicList<label> removedCells_;

public:

    // Constructor
    changeMap()
    :
        index_(-1),
        pIndex_(-1),
        type_(-1),
        addedPoints_(5),
        addedEdges_(5),
        addedFaces_(5),
        addedCells_(5),
        removedPoints_(5),
        removedEdges_(5),
        removedFaces_(5),
        removedCells_(5)
    {}

    //- Access

        // Entity index
        inline label& index();
        inline label index() const;

        // Coupled patch index
        inline label& patchIndex();
        inline label patchIndex() const;

        // Type
        inline label& type();
        inline label type() const;

        // Added entities
        inline void addPoint
        (
            const label pIndex,
            const labelList& master = labelList()
        );

        inline void addEdge
        (
            const label eIndex,
            const labelList& master = labelList()
        );

        inline void addFace
        (
            const label fIndex,
            const labelList& master = labelList()
        );

        inline void addCell
        (
            const label cIndex,
            const labelList& master = labelList()
        );

        // Return the list of added entities
        inline const List<objectMap>& addedPointList() const;
        inline const List<objectMap>& addedEdgeList() const;
        inline const List<objectMap>& addedFaceList() const;
        inline const List<objectMap>& addedCellList() const;

        // Removed entities
        inline void removePoint(const label pIndex);
        inline void removeEdge(const label eIndex);
        inline void removeFace(const label fIndex);
        inline void removeCell(const label cIndex);

        // Return the list of removed entities
        inline const labelList& removedPointList() const;
        inline const labelList& removedEdgeList() const;
        inline const labelList& removedFaceList() const;
        inline const labelList& removedCellList() const;

    //- Edit

        // Clear existing lists
        inline void clear();

    //- Operators

        inline void operator=(const changeMap& rhs);

    //- IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const changeMap&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "changeMapI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
