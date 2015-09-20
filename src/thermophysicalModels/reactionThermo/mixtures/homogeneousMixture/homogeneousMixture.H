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
    Foam::homogeneousMixture

Description
    Foam::homogeneousMixture

SourceFiles
    homogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousMixture_H
#define homogeneousMixture_H

#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class homogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class homogeneousMixture
:
    public basicMultiComponentMixture
{
    // Private data

        static const int nSpecies_ = 1;
        static const char* specieNames_[1];

        ThermoType reactants_;
        ThermoType products_;

        mutable ThermoType mixture_;

        //- Construct as copy (not implemented)
        homogeneousMixture(const homogeneousMixture<ThermoType>&);

        //- Regress variable
        volScalarField& b_;


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        homogeneousMixture(const dictionary&, const fvMesh&, const objectRegistry&);


    //- Destructor
    virtual ~homogeneousMixture()
    {}


    // Member functions

        const ThermoType& mixture(const scalar) const;

        const ThermoType& cellMixture(const label celli) const
        {
            return mixture(b_[celli]);
        }

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture(b_.boundaryField()[patchi][facei]);
        }

        const ThermoType& cellReactants(const label) const
        {
            return reactants_;
        }

        const ThermoType& patchFaceReactants(const label, const label) const
        {
            return reactants_;
        }

        const ThermoType& cellProducts(const label) const
        {
            return products_;
        }

        const ThermoType& patchFaceProducts(const label, const label) const
        {
            return products_;
        }

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "homogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
