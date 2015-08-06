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
    VRWGraphList

Description
    This class is an implementation of a list of graphs
    with variable column width. The implementation is memory efficient.

SourceFiles
    VRWGraphListI.H
    VRWGraphList.C

\*---------------------------------------------------------------------------*/

#ifndef VRWGraphList_H
#define VRWGraphList_H

#include "VRWGraph.H"
#include "subGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class VRWGraphList Declaration
\*---------------------------------------------------------------------------*/

class VRWGraphList
{
        
    // Private data
        //- graph containing the data
        VRWGraph data_;
    
        //- number of rows
        LongList<rowElement> rows_;

    // Private member functions
        //- check index
        inline void checkIndex
        (
            const label i,
            const label j,
            const label k
        ) const;
    
public:

    // Constructors

        //- Construct null
        inline VRWGraphList();

        //- Copy contructor
        inline VRWGraphList(const VRWGraphList&);

    // Destructor

        inline ~VRWGraphList();

    // Member Functions

        // Access

            //- Returns the number of graphs
            inline label size() const;
            
            //- Returns the number of rows in the graph at that position
            inline label sizeOfGraph(const label posI) const;
            
            //- Return the number of element in the row at the given position
            inline label sizeOfRow(const label posI, const label rowI) const;

            //- Clear the graph
            inline void clear();

    // Member Operators

        //- Append a graph at the end of the graphList
        template<class GraphType>
        inline void appendGraph(const GraphType& l);
        
        //- get and set operators
        inline label operator()
        (
            const label i,
            const label j,
            const label k
        ) const;
        
        inline label& operator()(const label i, const label j, const label k);
        
        inline const subGraph<const VRWGraph> operator[](const label i) const;
        
        //- Assignment operator
        inline void operator=(const VRWGraphList&);


    // IOstream operators

        // Write VRWGraphList to Ostream.
        friend Ostream& operator<<(Ostream&, const VRWGraphList&);

        //- Read from Istream, discarding contents of existing VRWGraphList.
/*        friend Istream& operator>> <T, width>
        (
            Istream&,
            VRWGraphList<T, width>&
        );
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VRWGraphListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
