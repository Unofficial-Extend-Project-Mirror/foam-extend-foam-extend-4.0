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
    Foam::dynamicMeshLoadBalance

Description
    Dynamic mesh load balance class.
    The modifier will examine the weighted imbalance in the mesh to determine
    when the load distribution is no longer satisfactory.
    Imbalance is measured as a weighted field over all cells and is checked over
    a parallel decomposition
    When load balance is detected:
    - parallel-aware decompositionMethod to receive the cell-to-processor
    assignment on each processor
    - polyMesh is split into parts and point, face, cell, boundary addressing
    objects are created
    - tools equivalent to decomposePar are used to split up the original mesh
    into processor parts
    - for each processor piece:
        - the new mesh is built, ready to be transferred to a processor
        - fields are subset onto each processor mesh
        - mesh and fields are passed onto the target processor
    - ALL DECOMPOSITION PIECES ARE CLEARED AFTER SEND, AND RE-USED ONTO RECEIVE
    - processor pieces are received from all components
    - new mesh is built by inserting processor pieces in order of processors
    - for each processor piece:
        - if different from myProcNo, point, face, cell, boundary addressing
        objects are created
        - if equal to myProcNo, a mapPolyMesh object is assembled
    - once the complete mesh is built, updateMesh is called on the local
    processor.  This will re-map local data and re-size the fields
    - when completed, tools equivalent to reconstructPar are used to rebuild
    the parts of fields from other processors

SourceFiles
    dynamicMeshLoadBalance.C
    dynamicMeshLoadBalanceIO.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMeshLoadBalance_H
#define dynamicMeshLoadBalance_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                     Class dynamicMeshLoadBalance Declaration
\*---------------------------------------------------------------------------*/

class dynamicMeshLoadBalance
{
    // Private data

        //- Dynamic mesh reference
        dynamicFvMesh& mesh_;

        //- Decomposition dictionary
        IOdictionary dict_;

        //- Imbalance trigger 0 < imbalance < 1
        scalar imbalanceTrigger_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicMeshLoadBalance(const dynamicMeshLoadBalance&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicMeshLoadBalance&);


        //- Send mesh to other processor
            static void sendMesh
            (
                const fvMesh& mesh,
                Ostream& os
            );
    
public:

    // //- Runtime type information
    // TypeName("testLoadBalanceFvMesh");


    // Constructors

        //- Construct from dynamic mesh
        dynamicMeshLoadBalance(dynamicFvMesh&);


    //- Destructor
    ~dynamicMeshLoadBalance();


    // Member Functions

        //- Check balance.  Return true if re-balancing is required
        bool checkLoadBalance(const scalarField& weights) const;

        //- Update mesh
        autoPtr<mapPolyMesh> updateMesh() const;

        //- Update fields
        void updateFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
