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
    processorFaPatch

Description
    Processor patch.

SourceFiles
    processorFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorFaPatch_H
#define processorFaPatch_H

#include "coupledFaPatch.H"
#include "processorLduInterface.H"
// #include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class processorFaPatch Declaration
\*---------------------------------------------------------------------------*/

class processorFaPatch
:
    public coupledFaPatch,
    public processorLduInterface
{
    // Private data

        int myProcNo_;
        int neighbProcNo_;

        //- Processor-neighbbour patch edge centres
        vectorField neighbEdgeCentres_;

        //- Processor-neighbbour patch edge lengths
        vectorField neighbEdgeLengths_;

        //- Processor-neighbbour patch neighbour face centres
        vectorField neighbEdgeFaceCentres_;

        //- Corresponding neighbouring local point label for every local point
        //  (so localPoints()[i] == neighb.localPoints()[neighbPoints_[i]])
        mutable labelList* neighbPointsPtr_;

        //- The set of labels of the processor patch points which are
        //  non-global, i.e. present in this processor patch
        mutable labelList* nonGlobalPatchPointsPtr_;

protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(scalarField&) const;

        //- Find non-globa patch points
        void makeNonGlobalPatchPoints() const;

protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        void initGeometry();

        //- Calculate the patch geometry
        void calcGeometry();

        //- Initialise the patches for moving points
        void initMovePoints(const pointField&);

        //- Correct patches after moving points
        void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();

public:

    //- Runtime type information
//     TypeName(processorPolyPatch::typeName_());
    TypeName("processor");


    // Constructors

        //- Construct from components
        processorFaPatch
        (
            const word& name,
            const labelList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label ngbPolyPatchIndex,
            const label myProcNo,
            const label neighbProcNo
        )
        :
            coupledFaPatch(name, edgeLabels, index, bm, ngbPolyPatchIndex),
            myProcNo_(myProcNo),
            neighbProcNo_(neighbProcNo),
            neighbEdgeCentres_(),
            neighbEdgeLengths_(),
            neighbEdgeFaceCentres_(),
            neighbPointsPtr_(NULL),
            nonGlobalPatchPointsPtr_(NULL)
        {}

        //- Construct from dictionary
        processorFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            coupledFaPatch(name, dict, index, bm),
            myProcNo_(readLabel(dict.lookup("myProcNo"))),
            neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
            neighbEdgeCentres_(),
            neighbEdgeLengths_(),
            neighbEdgeFaceCentres_(),
            neighbPointsPtr_(NULL),
            nonGlobalPatchPointsPtr_(NULL)
        {}

    // Destructor

        virtual ~processorFaPatch();


    // Member functions

        //- Return processor number
        int myProcNo() const
        {
            return myProcNo_;
        }

        //- Return neigbour processor number
        int neighbProcNo() const
        {
            return neighbProcNo_;
        }

        //- Return true if running parallel
        virtual bool coupled() const
        {
            if (Pstream::parRun())
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const
        {
            return coupledFaPatch::forwardT();
        }

        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const;


        //- Return processor-neighbbour patch edge centres
        const vectorField& neighbEdgeCentres() const
        {
            return neighbEdgeCentres_;
        }

        //- Return processor-neighbbour patch edge lengths
        const vectorField& neighbEdgeLengths() const
        {
            return neighbEdgeLengths_;
        }

        //- Return processor-neighbbour patch neighbour face centres
        const vectorField& neighbEdgeFaceCentres() const
        {
            return neighbEdgeFaceCentres_;
        }

        //- Return neighbour point labels. This is for my local point the
        //  corresponding local point on the other side. Call
        //  faBoundaryMesh::updateMesh() on all processors
        //  before using this.
        const labelList& neighbPoints() const;

        //- Return the set of labels of the processor patch points which are
        //  non-global, i.e. present in this processorFaPatch
        const labelList& nonGlobalPatchPoints() const;
    

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

            //- Write the patch data as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
