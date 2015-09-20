/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
#include "tolerancesSwitch.H"

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
        static const debug::tolerancesSwitch matchTol_;

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
