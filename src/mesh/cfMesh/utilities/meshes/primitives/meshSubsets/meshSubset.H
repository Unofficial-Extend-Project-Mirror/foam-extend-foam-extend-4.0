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
    meshSubset

Description
    A class containing labels of elements

SourceFiles
    meshSubsetI.H

\*---------------------------------------------------------------------------*/

#ifndef meshSubset_H
#define meshSubset_H

#include "label.H"
#include "HashSet.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Ostream;
class VRWGraph;

/*---------------------------------------------------------------------------*\
                        class meshSubset Declaration
\*---------------------------------------------------------------------------*/

class meshSubset
{
    // Private data
        //- name of the given subset
        word name_;

        //- type of subset
        label type_;

        //- labels of elements
        HashSet<label> data_;

public:

    // Enumerators

        enum subsetType_
        {
            UNKNOWN = 0,
            CELLSUBSET = 1,
            FACESUBSET = 2,
            POINTSUBSET = 4,
            FEATUREEDGESUBSET = 8
        };

    // Constructors
        //- Null constructor
        inline meshSubset();

        //- Construct from name. Creates an empty subset
        inline meshSubset(const word& name, const subsetType_&);

        //- Construct from name and element labels
        template<class ListType>
        inline meshSubset
        (
            const word& name,
            const subsetType_&,
            const ListType& elements
        );

        //- Copy construct
        inline meshSubset(const meshSubset&);

        //- Construct from Istream
        inline meshSubset(Istream&);

    // Destructor
        ~meshSubset();


    // Member Functions

        // Access

            //- Return name
            inline const word& name() const;

            //- Type of subset
            inline label type() const;

            //- elements contained in the subset
            template<class ListType>
            inline void containedElements(ListType&) const;

        // Edit

            //- add element label to subset
            inline void addElement(const label);

            //- remove element from subset
            inline void removeElement(const label);

            //- update subset after mesh modification
            //- the list contains new labels of mesh elements
            //- removed elements have a negative label
            template<class ListType>
            inline void updateSubset(const ListType&);

            //- update subset after modification
            //- this modifier is used in case when elements are decomposed
            inline void updateSubset(const VRWGraph&);

        // Searches
            //- find if the element exists in the subset
            inline bool contains(const label) const;

    // operators

            inline void operator=(const meshSubset&);
            inline bool operator==(const meshSubset&) const;
            inline bool operator!=(const meshSubset&) const;

    // IOstream operators

            friend inline Ostream& operator<<(Ostream&, const meshSubset&);
            friend inline Istream& operator>>(Istream&, meshSubset&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshSubsetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
