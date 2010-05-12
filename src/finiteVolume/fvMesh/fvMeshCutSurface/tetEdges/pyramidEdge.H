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
    Foam::pyramidEdge

Description
    Implicit description of pyramid edge (from tet decomposition).
    Pyramid edge is edge between mesh vertex and cell centre.

    Stores
        - vertex label
        - cell label

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef pyramidEdge_H
#define pyramidEdge_H

#include "label.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class pyramidEdge Declaration
\*---------------------------------------------------------------------------*/

class pyramidEdge
{
    // Private data

        //- vertex label
        label vertexLabel_;

        //- cell label
        label cellLabel_;

public:

    // Public classes

        //- Hash function
        class pyramidEdgeHash
        {

        public:

            pyramidEdgeHash()
            {}

            //- simple hashing function of labels
            label operator()(const pyramidEdge& pe, const label tableSize) const
            {
                // Note: mag since multiply might overflow and produce
                // negative numbers
                return
                    mag
                    (
                        pe.vertexLabel()
                      + (pe.vertexLabel() * pe.vertexLabel())
                      + pe.cellLabel()
                      + (pe.cellLabel() * pe.cellLabel())
                    )
                  % tableSize;
            }
        };


    // Constructors

        //- Construct null
        inline pyramidEdge()
        :
            vertexLabel_(-1),
            cellLabel_(-1)
        {}

        //- Construct from components
        inline pyramidEdge
        (
            const label vertexLabel,
            const label cellLabel
        )
        :
            vertexLabel_(vertexLabel),
            cellLabel_(cellLabel)
        {}

    // Member Functions

        label vertexLabel() const
        {   
            return vertexLabel_;
        }

        label cellLabel() const
        {   
            return cellLabel_;
        }

        template <class T>
        T interpolate
        (
            const Field<T>& cellField,
            const Field<T>& vertField,
            const scalar weight
        ) const
        {
            return
                (1-weight)*vertField[vertexLabel_]
              + weight*cellField[cellLabel_];
        }


        point coord(const primitiveMesh& mesh, const scalar weight) const
        {
            return interpolate(mesh.cellCentres(), mesh.points(), weight);
        }


    // Member Operators

        bool operator==(const pyramidEdge& pe) const
        {
            return
                (vertexLabel() == pe.vertexLabel())
             && (cellLabel() == pe.cellLabel());
        }


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream& os, const pyramidEdge& pe)
        {
            os  << token::BEGIN_LIST
                << pe.vertexLabel_ << token::SPACE
                << pe.cellLabel_
                << token::END_LIST;

            // Check state of Ostream
            os.check("Ostream& operator<<(Ostream&, const pyramidEdge&)");

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
