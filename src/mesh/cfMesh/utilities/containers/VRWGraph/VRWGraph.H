/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    VRWGraph

Description
    This class is an implementation of a graph with variable column width.
    The imoplementation is memory efficient.

SourceFiles
    VRWGraphI.H
    VRWGraph.C

\*---------------------------------------------------------------------------*/

#ifndef VRWGraph_H
#define VRWGraph_H

#include "labelLongList.H"
#include "graphRow.H"
#include "DynList.H"
#include "bool.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class VRWGraphModifier;

class rowElement
{
    // Private data
        //- starting element of the row
        label start_;

        //- number of elements in the row
        label size_;

    public:

        inline rowElement()
        :
            start_(),
            size_()
        {}

        inline rowElement(const label i, const label j)
        :
            start_(i),
            size_(j)
        {}

        inline ~rowElement()
        {}

        inline label start() const
        {
            return start_;
        }
        inline label& start()
        {
            return start_;
        }

        inline label size() const
        {
            return size_;
        }
        inline label& size()
        {
            return size_;
        }
};

/*---------------------------------------------------------------------------*\
                           Class VRWGraph Declaration
\*---------------------------------------------------------------------------*/

class VRWGraph
{
    // Private data
        //- list containing data
        labelLongList data_;

        //- number of rows
        LongList<rowElement> rows_;

    // Private member functions
        //- check index
        inline void checkIndex(const label i, const label j) const;

    // Enumerators
        enum typeOfEntries
        {
            NONE = 0,
            INVALIDROW=-10,
            FREEENTRY=-11,
            FREESTART=-12
        };

public:

    // Friend classes

        friend class VRWGraphSMPModifier;

    // Constructors

        //- Construct null
        inline VRWGraph();

        //- Construct given number of rows
        explicit inline VRWGraph(const label size);

        //- Construct given number of rows and row size
        explicit inline VRWGraph
        (
            const label nRows,
            const label nColumnsInRow
        );

        //- Construct to given number of rows, row size and initialize
        explicit inline VRWGraph
        (
            const label nRows,
            const label nColumnsInRow,
            const label t
        );

        //- Copy contructor
        inline VRWGraph(const VRWGraph&);

    // Destructor

        inline ~VRWGraph();

    // Member Functions

        // Access

            //- Returns the number of rows
            inline label size() const;

            //- Returns the number of elements in the given row
            inline label sizeOfRow(const label rowI) const;

        // Edit

            //- Reset the number of rows
            inline void setSize(const label);

            //- Reset the number of rows. The second argument specifies
            //- the reserved column width
            inline void setSizeAndColumnWidth
            (
                const label newNumRows,
                const label rcWidth
            );

            //- Set the number of rows and the size of each row
            template<class ListType>
            inline void setSizeAndRowSize(const ListType&);

            //- Reset the size of the given row
            inline void setRowSize(const label rowI, const label newSize);

            //- Clear the graph
            inline void clear();

    // Member Operators

        //- Append a list as a row at the end of the graph
        template<class ListType>
        inline void appendList(const ListType& l);

        //- Append an element to the given row
        inline void append(const label rowI, const label);

        //- Append an element to the given row if it does not exist there
        inline void appendIfNotIn(const label rowI, const label);

        //- Set row with the list
        template<class ListType>
        inline void setRow(const label rowI, const ListType& l);

        //- merge graphs with the identical number of rows
        //- into a single one. Use for SMP parallelisation
        inline void mergeGraphs(const List<VRWGraph>& graphParts);

        //- set the graph to the reverse of the original graph.
        //- the rows of such graph store the rows which contain the elements
        //- of the original graph
        template<class GraphType>
        inline void reverseAddressing
        (
            const label nRows,
            const GraphType& origGraph
        );

        template<class GraphType>
        inline void reverseAddressing(const GraphType& origGraph);

        inline void reverseAddressing
        (
            const label nRows,
            const VRWGraph& origGraph
        );

        inline void reverseAddressing(const VRWGraph& origGraph);

        //- optimize memory usage
        // this should be used once the graph will not be resized any more
        void optimizeMemoryUsage();

        //- check if the element is in the given row (takes linear time)
        inline bool contains(const label rowI, const label e) const;
        inline label containsAtPosition(const label rowI, const label e) const;

        //- get and set operators
        inline label operator()(const label i, const label j) const;
        inline label& operator()(const label i, const label j);

        inline constRow operator[](const label i) const;
        inline row operator[](const label i);

        //- Assignment operator
        inline void operator=(const VRWGraph&);


    // IOstream operators

        // Write VRWGraph to Ostream.
        friend Ostream& operator<<(Ostream&, const VRWGraph&);

        //- Read from Istream, discarding contents of existing VRWGraph.
/*        friend Istream& operator>> <T, width>
        (
            Istream&,
            VRWGraph<T, width>&
        );
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VRWGraphI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forAllRow(graph, rowI, index) \
    for(Foam::label index=0;index<(graph).sizeOfRow(rowI);++index)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
