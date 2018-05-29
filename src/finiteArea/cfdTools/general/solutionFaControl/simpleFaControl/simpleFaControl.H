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
    Foam::simpleFaControl

Description
    SIMPLE control class to supply convergence information/checks for
    the SIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef simpleFaControl_H
#define simpleFaControl_H

#include "solutionFaControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class simpleFaControl Declaration
\*---------------------------------------------------------------------------*/

class simpleFaControl
:
    public solutionFaControl
{

protected:

    // Protected Data

        //- Initialised flag
        bool initialised_;


    // Protected Member Functions

        //- Read controls from faSolution dictionary
        void read();

        //- Return true if all convergence checks are satisfied
        bool criteriaSatisfied();

        //- Disallow default bitwise copy construct
        simpleFaControl(const simpleFaControl&);

        //- Disallow default bitwise assignment
        void operator=(const simpleFaControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("simpleFaControl");


    // Constructors

        //- Construct from mesh
        simpleFaControl(faMesh& mesh);


    //- Destructor
    virtual ~simpleFaControl();


    // Member Functions

        // Solution control

            //- Loop loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
