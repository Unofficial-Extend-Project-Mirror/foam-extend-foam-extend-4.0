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
    verticalValveBank

Description
    A list of vertical valves.


\*---------------------------------------------------------------------------*/

#ifndef verticalValveBank_H
#define verticalValveBank_H

#include "PtrList.H"
#include "engineVerticalValve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class verticalValveBank Declaration
\*---------------------------------------------------------------------------*/

class verticalValveBank
:
    public PtrList<engineVerticalValve>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        verticalValveBank(const verticalValveBank&);

        //- Disallow default bitwise assignment
        void operator=(const verticalValveBank&);


public:

    // Constructors

        //- Construct from Istream
        verticalValveBank
        (
            const polyMesh& mesh,
            Istream& is
        )
        {
            PtrList<entry> valveEntries(is);

            setSize(valveEntries.size());

            forAll (valveEntries, valveI)
            {
                set
                (
                    valveI,
                    new engineVerticalValve
                    (
                        valveEntries[valveI].keyword(),
                        mesh,
                        valveEntries[valveI].dict()
                    )
                );
            }
        }


    // Destructor - default


    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
