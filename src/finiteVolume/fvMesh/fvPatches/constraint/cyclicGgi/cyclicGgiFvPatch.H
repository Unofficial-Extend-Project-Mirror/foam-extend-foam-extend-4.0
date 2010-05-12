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
    cyclicGgiFvPatch

Description
    Non-matching cyclic Generalized grid interface patch.

Author:
    Martin Beaudoin, Hydro-Quebec, (2008)

Contributor:
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    cyclicGgiFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGgiFvPatch_H
#define cyclicGgiFvPatch_H

#include "ggiFvPatch.H"
#include "cyclicGgiPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cyclicGgiFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicGgiFvPatch
:
    public ggiFvPatch
{
    // Private Data

        //- Reference to polyPatch
        const cyclicGgiPolyPatch& cyclicGgiPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        virtual void makeDeltaCoeffs(scalarField&) const;


public:

    //- Runtime type information
    TypeName(cyclicGgiPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        cyclicGgiFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            ggiFvPatch(patch, bm),
            cyclicGgiPolyPatch_(refCast<const cyclicGgiPolyPatch>(patch))
        {}


    // Destructor

        virtual ~cyclicGgiFvPatch()
        {}


    // Member functions

        // Access

            //- Return shadow patch
            const cyclicGgiFvPatch& shadow() const;

            //- Interpolate face field
            template<class Type>
            tmp<Field<Type> > interpolate(const Field<Type>& pf) const
            {
                return cyclicGgiPolyPatch_.interpolate(pf);
            }

            template<class Type>
            tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const
            {
                return cyclicGgiPolyPatch_.interpolate(tpf);
            }

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
