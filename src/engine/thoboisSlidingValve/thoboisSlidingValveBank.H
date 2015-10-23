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
    verticalValveBank

Description
    A list of vertical valves.


\*---------------------------------------------------------------------------*/

#ifndef thoboisSlidingValveBank_H
#define thoboisSlidingValveBank_H

#include "PtrList.H"
#include "thoboisSlidingValve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class verticalValveBank Declaration
\*---------------------------------------------------------------------------*/

class thoboisSlidingValveBank
:
    public PtrList<thoboisSlidingValve>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thoboisSlidingValveBank(const thoboisSlidingValveBank&);

        //- Disallow default bitwise assignment
        void operator=(const thoboisSlidingValveBank&);


public:

    // Constructors

        //- Construct from Istream
        thoboisSlidingValveBank
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
                    new thoboisSlidingValve
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
