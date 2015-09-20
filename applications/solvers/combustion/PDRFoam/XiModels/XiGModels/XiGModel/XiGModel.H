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
    Foam::XiGModel

Description
    Base-class for all Xi generation models used by the b-Xi combustion model.
    See Technical Report SH/RE/01R for details on the PDR modelling. For details
    on the use of XiGModel see \link XiModel.H \endlink. The model available is
    \link instabilityG.H \endlink

SourceFiles
    XiGModel.C

\*---------------------------------------------------------------------------*/

#ifndef XiGModel_H
#define XiGModel_H

#include "IOdictionary.H"
#include "hhuCombustionThermo.H"
#include "RASModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiGModel Declaration
\*---------------------------------------------------------------------------*/

class XiGModel
{

protected:

    // Protected data

        dictionary XiGModelCoeffs_;

        const hhuCombustionThermo& thermo_;
        const compressible::RASModel& turbulence_;
        const volScalarField& Su_;


private:

    // Private Member Functions

        //- Disallow copy construct
        XiGModel(const XiGModel&);

        //- Disallow default bitwise assignment
        void operator=(const XiGModel&);


public:

    //- Runtime type information
    TypeName("XiGModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiGModel,
            dictionary,
            (
                const dictionary& XiGProperties,
                const hhuCombustionThermo& thermo,
                const compressible::RASModel& turbulence,
                const volScalarField& Su
            ),
            (
                XiGProperties,
                thermo,
                turbulence,
                Su
            )
        );


    // Selectors

        //- Return a reference to the selected XiG model
        static autoPtr<XiGModel> New
        (
            const dictionary& XiGProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Constructors

        //- Construct from components
        XiGModel
        (
            const dictionary& XiGProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Destructor

        virtual ~XiGModel();


    // Member Functions

        //- Return the flame-wrinking genration rate
        virtual tmp<volScalarField> G() const = 0;

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const
        {
            return turbulence_.muEff();
        }

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiGProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
