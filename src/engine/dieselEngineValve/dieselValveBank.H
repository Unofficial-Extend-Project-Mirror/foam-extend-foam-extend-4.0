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
    dieselValveBank

Description
    A list of valves.


\*---------------------------------------------------------------------------*/

#ifndef dieselValveBank_H
#define dieselValveBank_H

#include "PtrList.H"
#include "dieselEngineValve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dieselValveBank Declaration
\*---------------------------------------------------------------------------*/

class dieselValveBank
:
    public PtrList<dieselEngineValve>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dieselValveBank(const dieselValveBank&);

        //- Disallow default bitwise assignment
        void operator=(const dieselValveBank&);


public:

    // Constructors

        //- Construct from Istream
        dieselValveBank
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
                    new dieselEngineValve
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
