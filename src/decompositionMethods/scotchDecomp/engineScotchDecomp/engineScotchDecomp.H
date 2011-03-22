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
    Foam::engineScotchDecomp

Description
    Domain decomposition for internal combustion engine simulations with
    topological changes.

    The algorithm attempts to be as general as possible and operates in the
    following steps:
    - identify piston region and peform a decomposition biased towards the
    piston axis.  This aims to preserve load balancing when cell layers
    are added
    - identify static regions and decompose them separately

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    engineScotchDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef engineScotchDecomp_H
#define engineScotchDecomp_H

#include "scotchDecomp.H"
#include "Switch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class engineScotchDecomp Declaration
\*---------------------------------------------------------------------------*/

class engineScotchDecomp
:
    public scotchDecomp
{
    // Private data

        //- Dictionary
        dictionary dict_;

        //- Sliding patch pairs
        List<Pair<word> > slidingPatchPairs_;

        //- Use expanded area around sliding patches
        Switch expandSliding_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        engineScotchDecomp(const engineScotchDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const engineScotchDecomp&);


        //- Check and print error message
        static void check(const int, const char*);

        label decompose
        (
            const List<int>& adjncy,
            const List<int>& xadj,
            const scalarField& cWeights,
            List<int>& finalDecomp
        );


public:

    //- Runtime type information
    TypeName("engineScotch");


    // Constructors

        //- Construct given the decomposition dictionary and mesh
        engineScotchDecomp
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Destructor

        virtual ~engineScotchDecomp()
        {}


    // Member Functions

        virtual bool parallelAware() const
        {
            return false;
        }

        //- Decompose cells with weights
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Decompose cells with weights with explicitly provided connectivity
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
