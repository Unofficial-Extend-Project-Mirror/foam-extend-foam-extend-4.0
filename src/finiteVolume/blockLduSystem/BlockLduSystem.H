/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::BlockLduSystem

Description
    BlockLduSystem is a wrapper for BlockLduMatrix with source field. This is
    the return type of implicit div and grad operators needed for implicitly
    coupled systems (namely p-U coupling).

Author
    Vuko Vukcevic, FMENA Zagreb.

SourceFiles
    BlockLduSystem.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduSystem_H
#define BlockLduSystem_H

#include "blockVectorNMatrices.H"
#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class blockType, class sourceType>
class BlockLduSystem;

template<class blockType, class sourceType>
Ostream& operator<<(Ostream&, const BlockLduSystem<blockType, sourceType>&);


/*---------------------------------------------------------------------------*\
                        Class BlockLduSystem Declaration
\*---------------------------------------------------------------------------*/

template<class blockType, class sourceType>
class BlockLduSystem
:
    public BlockLduMatrix<blockType>
{
    // Private data

        //- Source term
        Field<sourceType> source_;


public:

    // Constructors

        //- Construct given addressing
        explicit BlockLduSystem(const lduMesh&);

        //- Construct given addressing and source field
        BlockLduSystem(const lduMesh&, const Field<sourceType>&);

        //- Construct from components
        BlockLduSystem
        (
            const BlockLduMatrix<blockType>&,
            const Field<sourceType>&
        );

        //- Construct as copy
        BlockLduSystem(const BlockLduSystem<blockType, sourceType>&);


    // Destructor

       virtual ~BlockLduSystem()
       {}


    // Member functions

        //- Access

            Field<sourceType>& source()
            {
                return source_;
            }

            const Field<sourceType>& source() const
            {
                return source_;
            }


    // Member operators

        void operator=(const BlockLduSystem<blockType, sourceType>&);

        void negate();

        void operator+=(const BlockLduSystem<blockType, sourceType>&);
        void operator-=(const BlockLduSystem<blockType, sourceType>&);

        void operator*=(const scalarField&);
        void operator*=(const scalar);


    // Ostream operator

        friend Ostream& operator<< <blockType, sourceType>
        (
            Ostream&,
            const BlockLduSystem<blockType, sourceType>&
        );
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockLduSystem.C"
#   include "BlockLduSystemOperations.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
