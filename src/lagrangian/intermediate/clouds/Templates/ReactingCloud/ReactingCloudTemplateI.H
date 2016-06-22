/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::ReactingCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::multiComponentMixture<typename ParcelType::thermoType>&
Foam::ReactingCloud<ParcelType>::mcCarrierThermo() const
{
    return mcCarrierThermo_;
}


template<class ParcelType>
inline Foam::multiComponentMixture<typename ParcelType::thermoType>&
Foam::ReactingCloud<ParcelType>::mcCarrierThermo()
{
    return mcCarrierThermo_;
}


template<class ParcelType>
inline const Foam::CompositionModel<Foam::ReactingCloud<ParcelType> >&
Foam::ReactingCloud<ParcelType>::composition() const
{
    return compositionModel_;
}


template<class ParcelType>
inline const Foam::PhaseChangeModel<Foam::ReactingCloud<ParcelType> >&
Foam::ReactingCloud<ParcelType>::phaseChange() const
{
    return phaseChangeModel_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloud<ParcelType>::rhoTrans(const label i)
{
    return rhoTrans_[i];
}


template<class ParcelType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<ParcelType>::rhoTrans()
{
    return rhoTrans_;
}


template<class ParcelType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<ParcelType>::Srho(const label i) const
{
    return rhoTrans_[i]/(this->db().time().deltaT()*this->mesh().V());
}


template<class ParcelType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<ParcelType>::Srho() const
{
    tmp<DimensionedField<scalar, volMesh> > trhoTrans
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", rhoTrans_[0].dimensions(), 0.0)
        )
    );

    scalarField& sourceField = trhoTrans().field();
    forAll (rhoTrans_, i)
    {
        sourceField += rhoTrans_[i];
    }

    return trhoTrans/(this->db().time().deltaT()*this->mesh().V());
}


// ************************************************************************* //
