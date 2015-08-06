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
    subGraph

Description
    This class provides access to a row of a graph

SourceFiles
    subGraphI.H

\*---------------------------------------------------------------------------*/

#ifndef subGraph_H
#define subGraph_H

#include "bool.H"
#include "Ostream.H"
#include "error.H"
#include "graphRow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class graphType> class subGraph;
template<class graphType>
Ostream& operator<<(Ostream&, const subGraph<graphType>&);

/*---------------------------------------------------------------------------*\
                           Class subGraph Declaration
\*---------------------------------------------------------------------------*/

template<class graphType>
class subGraph
{
    // Private data
        //- reference to the graph
        graphType& data_;
    
        //- starts at row
        const label start_;
    
        //- number of rows in the subGraph
        const label size_;

    // Private member functions
        //- check index
        inline void checkIndex(const label i) const;
    
public:

    // Constructors

        //- Construct from graph and start row and size
        inline subGraph(graphType&, const label, const label);

        //- Copy contructor
        inline subGraph(const subGraph<graphType>&);

    // Destructor

        inline ~subGraph();

    // Member Functions

        //- Returns the number of rows
        inline label size() const;

        //- Returns the size of a given row
        inline label sizeOfRow(const label rowI) const;

    // Member Operators
        
        //- Append an element to the given row
        inline void append(const label rowI, const label);
        
        //- Append an element to the given row if it does not exist there
        inline void appendIfNotIn(const label rowI, const label);
        
        //- check if the element is in the given row (takes linear time)
        inline bool contains(const label rowI, const label e) const;
        inline label containsAtPosition(const label rowI, const label e) const;
        
        //- set and get operators
        inline label operator()(const label i, const label j) const;
        inline label& operator()(const label i, const label j);
        inline const graphRow<const graphType> operator[](const label) const;
        inline graphRow<graphType> operator[](const label);

    // IOstream operators

        // Write subGraph to Ostream.
        friend Ostream& operator<< <graphType>
        (
            Ostream&,
            const subGraph<graphType>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "subGraphI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
