/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::Cloud

Description

SourceFiles
    CloudDistributeTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef CloudDistributeTemplate_H
#define CloudDistributeTemplate_H

#include "CloudTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of functions
template<class ParticleType>
class Cloud;

class polyMesh;
class cloud;

/*---------------------------------------------------------------------------*\
                           Class Cloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class CloudDistribute
:
    public cloudDistribute
{
    // Private data

        //- Reference to the cloud
        Cloud<ParticleType>& cloud_;

        //- Transfer buffer for ordering by procID (for sending) and receiving
        List<IDLList<ParticleType> > transferList_;

public:

    // Constructors

        //- Construct from components and fill the transfer list
        CloudDistribute
        (
            const labelList& cellToProc,
            const labelListList& procCellAddressing,
            const labelListList& procFaceAddressing,
            cloud& cloud
        );


    // Member Functions

        //- Send particles in the transfer list
        void send(Ostream& toProc, const label procIndex);

        //- Receive particles into the transfer list
        void receive(Istream& fromProc, const label procIndex);

        //- Rebuild the cloud by adding the particles in the transfer list
        void rebuild
        (
            const PtrList<labelIOList>& cellProcAddressing,
            const PtrList<labelIOList>& faceProcAddressing
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CloudDistributeTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
