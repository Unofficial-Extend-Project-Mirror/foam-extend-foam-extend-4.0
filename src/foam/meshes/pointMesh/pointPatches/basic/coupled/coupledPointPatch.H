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
    Foam::coupledPointPatch

Description
    coupled patch for post-processing.  Used as the base class for processor
    and cyclic pointPatches

SourceFiles
    coupledPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledPointPatch_H
#define coupledPointPatch_H

#include "coupledPolyPatch.H"
#include "lduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointBoundaryMesh;

/*---------------------------------------------------------------------------*\
                      Class coupledPointPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledPointPatch
:
    // Bug fix: coupled patch must be derived from lduInterface
    // HJ, 16/Jan/2009
    public lduInterface
{
    // Private Member Functions

        //- Disallow default construct as copy
        coupledPointPatch(const coupledPointPatch&);

        //- Disallow default assignment
        void operator=(const coupledPointPatch&);


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry() = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry() = 0;

        //- Correct patches after moving points
        virtual void movePoints() = 0;

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh() = 0;

        //- Update of the patch topology
        virtual void updateMesh() = 0;


public:

    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        coupledPointPatch(const pointBoundaryMesh& bm);


    // Destructor

        virtual ~coupledPointPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return true;
            }


        // Dummy functionality for new lduInteface base class
        // To be removed.  HJ, 16/Jan/2009

            //- Return faceCell addressing
            virtual const unallocLabelList& faceCells() const
            {
                return labelList::null();
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const
            {
                return tmp<labelField>(NULL);
            }

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const
            {}

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const
            {
                return tmp<labelField>(NULL);
            }

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const
            {}

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const
            {
                return tmp<labelField>(NULL);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
