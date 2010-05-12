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
    tetPolyMeshLduAddressingCellDecomp

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef tetPolyMeshLduAddressingCellDecomp_H
#define tetPolyMeshLduAddressingCellDecomp_H

#include "lduAddressing.H"
#include "tetPolyMeshCellDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class tetPolyMeshLduAddressingCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyMeshLduAddressingCellDecomp
:
    public lduAddressing
{
    // Private data

        //- Lower
        labelList lowerAddr_;

        //- Upper
        labelList upperAddr_;

        //- Patch addressing
        labelListList patchAddr_;

        //- Patch field evaluation schedule
        const lduSchedule& patchSchedule_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyMeshLduAddressingCellDecomp
        (
            const tetPolyMeshLduAddressingCellDecomp&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPolyMeshLduAddressingCellDecomp&);


public:

    // Constructors

        //- Construct from components
        tetPolyMeshLduAddressingCellDecomp(const tetPolyMeshCellDecomp& mesh);


    // Destructor

        virtual ~tetPolyMeshLduAddressingCellDecomp()
        {}


    // Member Functions

        //- Return number of interfaces
        virtual label nPatches() const
        {
            return patchAddr_.size();
        }

        //- Return lower addressing (i.e. lower label = upper triangle)
        virtual const unallocLabelList& lowerAddr() const
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
        virtual const unallocLabelList& upperAddr() const
        {
            return upperAddr_;
        }

        //- Return patch addressing
        virtual const unallocLabelList& patchAddr(const label i) const
        {
            return patchAddr_[i];
        }

        // Return patch field evaluation schedule
        virtual const lduSchedule& patchSchedule() const
        {
            return patchSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
