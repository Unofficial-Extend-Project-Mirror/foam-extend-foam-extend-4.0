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
    Foam::oversetFvPatch

Description
    Overset patch.

SourceFiles
    oversetFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvPatch_H
#define oversetFvPatch_H

#include "coupledFvPatch.H"
#include "oversetLduInterface.H"
#include "oversetPolyPatch.H"
#include "oversetFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class oversetMesh;

/*---------------------------------------------------------------------------*\
                        Class oversetFvPatch Declaration
\*---------------------------------------------------------------------------*/

class oversetFvPatch
:
    public coupledFvPatch,
    public oversetLduInterface
{
    // Private Data

        //- Reference to overset patch
        const oversetPolyPatch& overPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(scalarField&) const;


public:

    //- Runtime type information
    TypeName(oversetPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        oversetFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            overPolyPatch_(refCast<const oversetPolyPatch>(patch))
        {}


    // Destructor

        virtual ~oversetFvPatch()
        {}


    // Member functions

        // Access

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            //- Return access to overset mesh
            const oversetMesh& overset() const;

            //- Return true: coupled patch
            virtual bool coupled() const
            {
                return true;
            }

        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& internalData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
