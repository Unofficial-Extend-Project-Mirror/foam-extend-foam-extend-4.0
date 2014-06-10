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
    graphRow

Description
    This class provides access to a row of a graph

SourceFiles
    graphRowI.H

\*---------------------------------------------------------------------------*/

#ifndef graphRow_H
#define graphRow_H

#include "bool.H"
#include "Ostream.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
class VRWGraph;

template<class graphType> class graphRow;
template<class graphType>
Ostream& operator<<(Ostream&, const graphRow<graphType>&);

/*---------------------------------------------------------------------------*\
                           Class graphRow Declaration
\*---------------------------------------------------------------------------*/

template<class graphType>
class graphRow
{
    // Private data
        //- reference to the graph
        graphType& data_;
    
        //- row number
        const label rowI_;

    // Private member functions
        //- check index
        inline void checkIndex(const label i) const;
    
public:

    // Constructors

        //- Construct from graph and row number
        inline graphRow(graphType&, const label);

        //- Copy contructor
        inline graphRow(const graphRow<graphType>&);

    // Destructor

        inline ~graphRow();

    // Member Functions

        //- Returns the number of rows
        inline label size() const;

        //- Reset the number of rows
        inline void setSize(const label size);

        //- Clear the graph
        inline void clear();

    // Member Operators
        
        //- Append an element to the given row
        inline void append(const label);
        
        //- Append an element to the given row if it does not exist there
        inline void appendIfNotIn(const label);
        
        //- check if the element is in the given row (takes linear time)
        inline bool contains(const label e) const;
        inline label containsAtPosition(const label e) const;
        
        //- set and get operators
        inline label operator[](const label) const;
        inline label& operator[](const label);
        
        //- Assignment operator
        inline void operator=(const graphRow<graphType>&);
        template<class listType>
        inline void operator=(const listType&);

    // IOstream operators

        // Write graphRow to Ostream.
        friend Ostream& operator<< <graphType>
        (
            Ostream&,
            const graphRow<graphType>&
        );

        //- Read from Istream, discarding contents of existing graphRow.
/*        friend Istream& operator>> <T, width>
        (
            Istream&,
            graphRow<T, width>&
        );
*/
};

typedef const graphRow<const VRWGraph> constRow;
typedef graphRow<VRWGraph> row;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "graphRowI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
