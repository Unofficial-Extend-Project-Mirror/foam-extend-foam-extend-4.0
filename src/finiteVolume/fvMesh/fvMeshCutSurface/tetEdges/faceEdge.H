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
    Foam::faceEdge

Description
    Implicit description of face edge (from tet decomposition).
    Face edge is edge between vertex of face and face centre.

    Stores
        - face label
        - index in face (gives the mesh vertex)


SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef faceEdge_H
#define faceEdge_H

#include "label.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class faceEdge Declaration
\*---------------------------------------------------------------------------*/

class faceEdge
{
    // Private data

        //- face label
        label faceLabel_;

        //- index in face
        label faceIndex_;

public:

    // Public classes

        //- Hash function
        class faceEdgeHash
        {

        public:

            faceEdgeHash()
            {}

            //- simple hashing function of labels
            label operator()(const faceEdge& fe, const label tableSize)
                const
            {
                // Note: mag since multiply might overflow and produce
                // negative numbers
                return
                    mag
                    (
                        fe.faceLabel()
                      + (fe.faceLabel() * fe.faceLabel())
                      + fe.faceIndex()
                      + (fe.faceIndex() * fe.faceIndex())
                    )
                  % tableSize;
            }
        };


    // Constructors

        //- Construct null
        inline faceEdge()
        :
            faceLabel_(-1),
            faceIndex_(-1)
        {}

        //- Construct from components
        inline faceEdge(const label faceLabel, const label faceIndex)
        :
            faceLabel_(faceLabel),
            faceIndex_(faceIndex)
        {}

    // Member Functions

        label faceLabel() const
        {
            return faceLabel_;
        }

        label faceIndex() const
        {
            return faceIndex_;
        }

        template <class T>
        T interpolate
        (
            const primitiveMesh& mesh,
            const Field<T>& faceField,
            const Field<T>& vertField,
            const scalar weight
        ) const
        {
            const face& f = mesh.faces()[faceLabel_];

            return
                (1-weight)*vertField[f[faceIndex_]]
              + weight*faceField[faceLabel_];
        }


        point coord(const primitiveMesh& mesh, const scalar weight) const
        {
            return interpolate(mesh, mesh.faceCentres(), mesh.points(), weight);
        }


    // Member Operators

        bool operator==(const faceEdge& fe) const
        {
            return
                (faceLabel() == fe.faceLabel())
             && (faceIndex() == fe.faceIndex());
        }


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream& os, const faceEdge& fe)
        {
            os  << token::BEGIN_LIST
                << fe.faceLabel_ << token::SPACE
                << fe.faceIndex_
                << token::END_LIST;

            // Check state of Ostream
            os.check("Ostream& operator<<(Ostream&, const faceEdge&)");

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
