/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    yamadaMises

Description
    Stress return method similar to that described by a number of
    references:
    Yamada, Y. et al. 1968. with Yoshimura, N. and Sakurai, T. Plastic
    stress-strain matrix and its application for the solution of
    elastic-plastic problems by the finite element method.
    International Journal of Mechanical Science. vol. 10, p. 343-354.
    and
    Maneeratana, K. 2000. Development of the finite volume method for
    non-linear structural applications.

    This approach gives the stress rate explicitly as a function of
    strain rates and instantaneous stresses

SourceFiles
    yamadaMises.C

Author
    Aleksandar Karac Zenica/UCD
    Rearranged by Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef yamadaMises_H
#define yamadaMises_H

#include "plasticityStressReturn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  class constitutiveModel;

/*---------------------------------------------------------------------------*\
                         Class yamadaMises Declaration
\*---------------------------------------------------------------------------*/

class yamadaMises
:
    public plasticityStressReturn
{
    // Private data

        // Reference to constitutive model
        constitutiveModel& constitutiveModel_;

        // Yield stress
        volScalarField sigmaY_;

        // Increment of plastic strain
        volSymmTensorField DEpsilonP_;

        // Increment of equivalent plastic strain
        volScalarField DEpsilonPEq_;

        // Total equivalent plastic strain
        volScalarField epsilonPEq_;

        // Active yielding flag
        // 1.0 for active yielding
        // 0.0 for fully elastic
        // 0.0 < beta < 1.0 for step which are a combination of fully elastic
        // and elastic-plastic
        volScalarField beta_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        yamadaMises(const yamadaMises&);

        //- Disallow default bitwise assignment
        void operator=(const yamadaMises&);

public:

    //- Runtime type information
    TypeName("yamadaMises");

    // Static data members


    // Constructors

        //- Construct from dictionary
        yamadaMises
        (
	 const word& name,
	 constitutiveModel& constitutiveModel
        );


    // Destructor

        virtual ~yamadaMises();


    // Member Functions

        //- Correct the rheological model
        virtual void correct();

        //- Return slave friction value fraction
        virtual const volSymmTensorField& DEpsilonP() const
        {
	  return DEpsilonP_;
	}

        //- Return slave friction value fraction
        virtual void updateYieldStress();

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
