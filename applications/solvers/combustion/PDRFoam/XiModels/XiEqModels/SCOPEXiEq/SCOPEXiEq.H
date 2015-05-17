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
    Foam::XiEqModel::SCOPEXiEq

Description
    Simple SCOPEXiEq model for XiEq based on SCOPEXiEqs correlation
    with a linear correction function to give a plausible profile for XiEq.
    See \link SCOPELaminarFlameSpeed.H \endlink for details on the SCOPE laminar
    flame speed model.

SourceFiles
    SCOPEXiEq.C

\*---------------------------------------------------------------------------*/

#ifndef SCOPEXiEq_H
#define SCOPEXiEq_H

#include "XiEqModel.H"
#include "SCOPELaminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class SCOPEXiEq Declaration
\*---------------------------------------------------------------------------*/

class SCOPEXiEq
:
    public XiEqModel
{
    // Private data

        scalar XiEqCoef;
        scalar XiEqExp;
        scalar lCoef;
        dimensionedScalar SuMin;

        //- The SCOPE laminar flame speed model used to obtain the
        //  Marstein number.  Note: the laminar flame speed need not be
        //  obtained form the same model.
        laminarFlameSpeedModels::SCOPE MaModel;


    // Private Member Functions

        //- Disallow copy construct
        SCOPEXiEq(const SCOPEXiEq&);

        //- Disallow default bitwise assignment
        void operator=(const SCOPEXiEq&);


public:

    //- Runtime type information
    TypeName("SCOPEXiEq");


    // Constructors

        //- Construct from components
        SCOPEXiEq
        (
            const dictionary& XiEqProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Destructor

        virtual ~SCOPEXiEq();


    // Member Functions

        //- Return the flame-wrinking XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
