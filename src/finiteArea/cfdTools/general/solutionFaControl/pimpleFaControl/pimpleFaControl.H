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
    Foam::pimpleFaControl

Description
    PIMPLE control class to supply convergence information/checks for
    the PIMPLE loop.

    May also be used to for PISO-based algorithms as PISO controls are a
    sub-set of PIMPLE controls.

\*---------------------------------------------------------------------------*/

#ifndef pimpleFaControl_H
#define pimpleFaControl_H

#include "solutionFaControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pimpleFaControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleFaControl
:
    public solutionFaControl
{
    // Private member functions

        //- Disallow default bitwise copy construct
        pimpleFaControl(const pimpleFaControl&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleFaControl&);


protected:

    // Protected data

        // Solution controls

            //- Maximum number of PIMPLE correctors
            label nCorrPIMPLE_;

            //- Maximum number of PISO correctors
            label nCorrPISO_;

            //- Current PISO corrector
            label corrPISO_;

            //- Flag to indicate whether to only solve turbulence on final iter
            bool turbOnFinalIterOnly_;

            //- Converged flag
            bool converged_;


    // Protected Member Functions

        //- Read controls from faSolution dictionary
        virtual void read();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pimpleFaControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        pimpleFaControl(faMesh& mesh, const word& dictName="PIMPLE");


    //- Destructor
    virtual ~pimpleFaControl();


    // Member Functions

        // Access

            //- Maximum number of PIMPLE correctors
            inline label nCorrPIMPLE() const;

            //- Maximum number of PISO correctors
            inline label nCorrPISO() const;

            //- Current PISO corrector index
            inline label corrPISO() const;


        // Solution control

            //- PIMPLE loop
            virtual bool loop();

            //- Pressure corrector loop
            inline bool correct();

            //- Helper function to identify when to store the intial residuals
            inline bool storeInitialResiduals() const;

            //- Helper function to identify first PIMPLE (outer) iteration
            inline bool firstIter() const;

            //- Helper function to identify final PIMPLE (outer) iteration
            inline bool finalIter() const;

            //- Helper function to identify final inner iteration
            inline bool finalInnerIter() const;

            //- Helper function to identify whether to solve for turbulence
            inline bool turbCorr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleFaControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
