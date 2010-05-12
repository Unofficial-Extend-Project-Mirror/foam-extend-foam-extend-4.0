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
    Foam::diagonalEdge

Description
    Implicit description of diagonal edge (from tet decomposition).
    Diangonal edge is edge between two vertices of same face.

    Stores
        - face label
        - indices in face

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef diagonalEdge_H
#define diagonalEdge_H

#include "label.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class diagonalEdge Declaration
\*---------------------------------------------------------------------------*/

class diagonalEdge
{
    // Private data

        //- face label
        label faceLabel_;

        //- index in face
        label faceIndex0_;

        //- index in face
        label faceIndex1_;

public:

    // Public classes

        //- Hash function
        class diagonalEdgeHash
        {

        public:

            diagonalEdgeHash()
            {}

            //- simple hashing function of labels. Commutative on face indices
            //  since edge is same.
            label operator()(const diagonalEdge& de, const label tableSize)
                const
            {
                // Note: mag since multiply might overflow and produce
                // negative numbers
                return
                    mag
                    (
                        de.faceLabel()
                      + (de.faceLabel()*de.faceLabel())
                      + de.faceIndex0() * de.faceIndex1()
                      + (de.faceIndex0() + de.faceIndex1())
                    )
                  % tableSize;
            }
        };


    // Constructors

        //- Construct null
        inline diagonalEdge()
        :
            faceLabel_(-1),
            faceIndex0_(-1),
            faceIndex1_(-1)
        {}

        //- Construct from components
        inline diagonalEdge
        (
            const label faceLabel,
            const label faceIndex0,
            const label faceIndex1
        )
        :
            faceLabel_(faceLabel),
            faceIndex0_(faceIndex0),
            faceIndex1_(faceIndex1)
        {}

    // Member Functions

        label faceLabel() const
        {   
            return faceLabel_;
        }

        label faceIndex0() const
        {   
            return faceIndex0_;
        }

        label faceIndex1() const
        {   
            return faceIndex1_;
        }

        template <class T>
        T interpolate
        (
            const primitiveMesh& mesh,
            const Field<T>& vertField,
            const scalar weight
        ) const
        {
            const face& f = mesh.faces()[faceLabel_];

            return
                (1-weight)*vertField[f[faceIndex0_]]
              + weight*vertField[f[faceIndex1_]];
        }

        point coord(const primitiveMesh& mesh, const scalar weight) const
        {
            return interpolate(mesh, mesh.points(), weight);
        }


    // Member Operators

        bool operator==(const diagonalEdge& de) const
        {
            return
                (faceLabel() == de.faceLabel())
             && (faceIndex0() == de.faceIndex0())
             && (faceIndex1() == de.faceIndex1());
        }


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream& os, const diagonalEdge& de)
        {
            os  << token::BEGIN_LIST
                << de.faceLabel_ << token::SPACE
                << de.faceIndex0_ << token::SPACE
                << de.faceIndex1_
                << token::END_LIST;

            // Check state of Ostream
            os.check("Ostream& operator<<(Ostream&, const diagonalEdge&)");

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
