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
    Foam::XiEqModels::basicSubGrid

Description
    Basic sub-grid obstacle flame-wrinking enhancement factor model.
    Details supplied by J Puttock 2/7/06.

    <b> Sub-grid flame area generation </b>

    \f$ n = N - \hat{\dwea{\vec{U}}}.n_{s}.\hat{\dwea{\vec{U}}} \f$
    \f$ n_{r} = \sqrt{n} \f$

    where:

        \f$ \hat{\dwea{\vec{U}}} = \dwea{\vec{U}} / \vert \dwea{\vec{U}}
        \vert \f$

        \f$ b = \hat{\dwea{\vec{U}}}.B.\hat{\dwea{\vec{U}}} / n_{r} \f$

    where:

        \f$ B \f$ is the file "B".

        \f$ N \f$ is the file "N".

        \f$  n_{s} \f$ is the file "ns".

    The flame area enhancement factor \f$ \Xi_{sub} \f$ is expected to
    approach:

    \f[
        \Xi_{{sub}_{eq}} =
            1 + max(2.2 \sqrt{b}, min(0.34 \frac{\vert \dwea{\vec{U}}
            \vert}{{\vec{U}}^{'}}, 1.6)) \times min(\frac{n}{4}, 1)
    \f]


SourceFiles
    basicSubGrid.C

\*---------------------------------------------------------------------------*/

#ifndef basicSubGrid_H
#define basicSubGrid_H

#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class basicSubGrid Declaration
\*---------------------------------------------------------------------------*/

class basicSubGrid
:
    public XiEqModel
{
    // Private data

        //- Count
        volScalarField N_;

        //- Sub-count
        volSymmTensorField ns_;

        //- tblock
        volSymmTensorField B_;

        //- Typical obstacle diameters per cell
        volScalarField Lobs_;

        //- Equilibrium Xi model due to turbulence
        autoPtr<XiEqModel> XiEqModel_;


    // Private Member Functions

        //- Disallow copy construct
        basicSubGrid(const basicSubGrid&);

        //- Disallow default bitwise assignment
        void operator=(const basicSubGrid&);


public:

    //- Runtime type information
    TypeName("basicSubGrid");


    // Constructors

        //- Construct from components
        basicSubGrid
        (
            const dictionary& XiEqProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Destructor

        virtual ~basicSubGrid();


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
