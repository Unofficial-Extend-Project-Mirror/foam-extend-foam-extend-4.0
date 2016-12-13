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
    Foam::mapDistributeLagrangian

Description
    Class containing mesh-to-mesh mapping information for particles

SourceFiles
    mapDistributeLagrangian.C

\*---------------------------------------------------------------------------*/

#ifndef mapDistributeLagrangian_H
#define mapDistributeLagrangian_H

#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class mapDistributeLagrangian Declaration
\*---------------------------------------------------------------------------*/

class mapDistributeLagrangian
{
    // Private data

        //- Map to distribute particles
        const mapDistribute particleMap_;

        //- Per element in subsetted mesh the cell label
        const labelListList constructCellLabels_;


public:

    // Constructors

        //- Construct from components
        mapDistributeLagrangian
        (
            const label nNewParticles,
            const labelListList& subParticleMap,
            const labelListList& constructParticleMap,
            const labelListList& constructCellLabels
        )
        :
            particleMap_(nNewParticles, subParticleMap, constructParticleMap),
            constructCellLabels_(constructCellLabels)
        {}

        //- Construct from components and steal storage
        mapDistributeLagrangian
        (
            const label nNewParticles,
            labelListList& subParticleMap,
            labelListList& constructParticleMap,
            labelListList& constructCellLabels,
            const bool reUse
        )
        :
            particleMap_
            (
                nNewParticles,
                subParticleMap,
                constructParticleMap,
                reUse
            ),
            constructCellLabels_(constructCellLabels, reUse)
        {}


    // Member Functions

        // Access

            //- Distribution map
            const mapDistribute& particleMap() const
            {
                return particleMap_;
            }

            //- Per received particle the destination cell label
            const labelListList& constructCellLabels() const
            {
                return constructCellLabels_;
            }


        // Edit

            //- distribute list of lagrangian data
            template<class T>
            void distributeLagrangianData(List<T>& lst) const
            {
                particleMap_.distribute(lst);
            }

            //- Correct for topo change.
            void updateMesh(const mapPolyMesh&)
            {
                notImplemented
                (
                    "mapDistributeLagrangian::updateMesh(const mapPolyMesh&)"
                );
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
