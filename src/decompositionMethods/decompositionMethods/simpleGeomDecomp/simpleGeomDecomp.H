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
    Foam::simpleGeomDecomp

Description

SourceFiles
    simpleGeomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef simpleGeomDecomp_H
#define simpleGeomDecomp_H

#include "geomDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class simpleGeomDecomp Declaration
\*---------------------------------------------------------------------------*/

class simpleGeomDecomp
:
    public geomDecomp
{
    // Private Member Functions

        void assignToProcessorGroup(labelList& processorGroup, const label);

        void assignToProcessorGroup
        (
            labelList& processorGroup,
            const label nProcGroup,
            const labelList& indices,
            const scalarField& weights,
            const scalar summedWeights
        );

        //- Disallow default bitwise copy construct
        simpleGeomDecomp(const simpleGeomDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const simpleGeomDecomp&);


public:

    //- Runtime type information
    TypeName("simple");


    // Constructors

        //- Construct given the decomposition dictionary
        explicit simpleGeomDecomp(const dictionary& decompositionDict);

        //- Construct given the decomposition dictionary and mesh
        simpleGeomDecomp
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Destructor

        virtual ~simpleGeomDecomp()
        {}


    // Member Functions

        virtual bool parallelAware() const
        {
            // simpleDecomp does not attempt to do anything across proc
            // boundaries
            return false;
        }

        virtual labelList decompose
        (
            const pointField& points
        );

        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        )
        {
            return decompose(cc, cWeights);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
