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
    Foam::cyclicFaPatch

Description
    Cyclic-plane patch.

SourceFiles
    cyclicFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFaPatch_H
#define cyclicFaPatch_H

#include "coupledFaPatch.H"
#include "cyclicLduInterface.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFaPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicFaPatch
:
    public coupledFaPatch,
    public cyclicLduInterface
{
    // Private data

    // Private member functions

    void calcTransforms();

protected:

    // Protected static data

        //- Relative tolerance (for geometric matching). Is factor of
        //  maximum edge length per face.
        static const scalar matchTol_;

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(scalarField&) const;


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from dictionary
        cyclicFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            coupledFaPatch(name, dict, index, bm)
        {}


    // Destructor

        virtual ~cyclicFaPatch()
        {}


    // Member functions

        // Access

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return coupledFaPatch::forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return coupledFaPatch::reverseT();
            }

            //- Initialise the calculation of the patch geometry
            virtual void initGeometry();

            //- Calculate the patch geometry
            virtual void calcGeometry();

            //- Initialise the patches for moving points
            virtual void initMovePoints(const pointField&);

            //- Correct patches after moving points
            virtual void movePoints(const pointField&);

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
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
