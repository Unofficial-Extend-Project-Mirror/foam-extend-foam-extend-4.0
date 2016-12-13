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

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#include "realGasHThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MixtureType>
void Foam::realGasHThermo<MixtureType>::calculate()
{
    const scalarField& hCells = h_.internalField();
    const scalarField& pCells = this->p_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& rhoCells= this->rho_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& drhodhCells = this->drhodh_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        mixture_.TH(hCells[celli], TCells[celli], pCells[celli], rhoCells[celli]);
        psiCells[celli] = mixture_.psi(rhoCells[celli], TCells[celli]);
        drhodhCells[celli] = mixture_.drhodH(rhoCells[celli], TCells[celli]);
        muCells[celli] = mixture_.mu(TCells[celli]);
        alphaCells[celli] = mixture_.alpha(rhoCells[celli], TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];
        fvPatchScalarField& pdrhodh = this->drhodh_.boundaryField()[patchi];
        fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        fvPatchScalarField& ph = h_.boundaryField()[patchi];
        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                prho[facei] = mixture_.rho(pp[facei], pT[facei],prho[facei]);
                ppsi[facei]=mixture_.psi(prho[facei],pT[facei]);
                pdrhodh[facei]=mixture_.drhodH(prho[facei],pT[facei]);
                ph[facei] = mixture_.H(prho[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pT[facei]);
                palpha[facei] = mixture_.alpha(prho[facei],pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                mixture_.TH(ph[facei], pT[facei],pp[facei],prho[facei]);
                pmu[facei] = mixture_.mu(pT[facei]);
                ppsi[facei]=mixture_.psi(prho[facei],pT[facei]);
                pdrhodh[facei]=mixture_.drhodH(prho[facei],pT[facei]);
                palpha[facei] = mixture_.alpha(prho[facei],pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::realGasHThermo<MixtureType>::realGasHThermo
(
    const fvMesh& mesh,
    const objectRegistry& obj
)
:
    basicPsiThermo(mesh, obj),
    MixtureType(*this, mesh, obj),

    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            obj,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->hBoundaryTypes()
    ),
    rho_
    (
        IOobject
        (
            "rhoThermo",
            mesh.time().timeName(),
            obj,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),
    drhodh_
    (
        IOobject
        (
            "drhodh",
            mesh.time().timeName(),
            obj,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -5, 2, 0, 0)
    )
{
    scalarField& hCells = h_.internalField();
    const scalarField& TCells = this->T_.internalField();
    const scalarField& pCells =this->p_.internalField();
    scalarField& rhoCells =this->rho_.internalField();

    forAll(rhoCells, celli)
    {
        rhoCells[celli]=this->cellMixture(celli).rho(pCells[celli],TCells[celli]);
    }

    forAll(rho_.boundaryField(), patchi)
    {
        rho_.boundaryField()[patchi] ==
            rho(this->T_.boundaryField()[patchi], patchi);
    }

    forAll(hCells, celli)
    {
        hCells[celli] = this->cellMixture(celli).H(rhoCells[celli],TCells[celli]);
    }

    forAll(h_.boundaryField(), patchi)
    {
        h_.boundaryField()[patchi] ==
            h(this->T_.boundaryField()[patchi], patchi);
    }

    hBoundaryCorrection(h_);

    calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::realGasHThermo<MixtureType>::~realGasHThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::realGasHThermo<MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering realGasHThermo<MixtureType>::correct()" << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "exiting realGasHThermo<MixtureType>::correct()" << endl;
    }
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::realGasHThermo<MixtureType>::h
(
    const scalarField& T,
    const labelList& cells
) const
{
    //CL: need the pressure of the internal field to calculate the realGas enthalpy
    //CL: this is done this way to assure compatibility to old OF Thermo-Versions
    const scalarField& pCells = this->p_.internalField();

    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, celli)
    {
        h[celli] = this->cellMixture(cells[celli]).H(this->cellMixture(cells[celli]).rho(pCells[cells[celli]],T[celli]),T[celli]);
    }

    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::realGasHThermo<MixtureType>::h
(
    const scalarField& T,
    const label patchi
) const
{
    //CL: need the pressure at the patch to calculate the realGas enthalpy
    //CL: this is done this way to assure compatibility to old OF Thermo-Versions
    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];

    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, facei)
    {
        h[facei] = this->patchFaceMixture(patchi, facei).H(this->patchFaceMixture(patchi, facei).rho(pp[facei], T[facei]),T[facei]);
    }

    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::realGasHThermo<MixtureType>::rho
(
    const scalarField& T,
    const label patchi
) const
{
    //CL: need the pressure at the patch to calculate the realGas enthalpy
    //CL: this is done this way to assure compatibility to old OF Thermo-Versions
    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];

    tmp<scalarField> trho(new scalarField(T.size()));
    scalarField& rho = trho();

    forAll(T, facei)
    {
        rho[facei] = this->patchFaceMixture(patchi, facei).rho(pp[facei], T[facei]);
    }

    return trho;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::realGasHThermo<MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    //CL: need the pressure at the patch to calculate the realGas enthalpy
    //CL: this is done this way to assure compatibility to old OF Thermo-Versions
    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];

    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    forAll(T, facei)
    {
        cp[facei] = this->patchFaceMixture(patchi, facei).Cp(this->patchFaceMixture(patchi, facei).rho(pp[facei], T[facei]),T[facei]);
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::realGasHThermo<MixtureType>::Cp() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                this->T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cp = tCp();

    forAll(this->T_, celli)
    {
        cp[celli] = this->cellMixture(celli).Cp(this->rho_[celli], this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cp.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pCp[facei] = this->patchFaceMixture(patchi, facei).Cp(prho[facei], pT[facei]);
        }
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::realGasHThermo<MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    //CL: need the pressure at the patch to calculate the realGas enthalpy
    //CL: this is done this way to assure compatibility to old OF Thermo-Versions
    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];

    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& cv = tCv();

    forAll(T, facei)
    {
        cv[facei] = this->patchFaceMixture(patchi, facei).Cv(this->patchFaceMixture(patchi, facei).rho(pp[facei], T[facei]), T[facei]);
    }

    return tCv;
}


// CL: Maybe this function should be changed so that it is not "const" fucntion anymore
template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::realGasHThermo<MixtureType>::rho()  const
{

    const fvMesh& mesh = this->T_.mesh();

    //CL: create a rho Field, which will be return
    //CL: the problem is that this function is "const",
    //CL: so a new variabel is needed
      tmp<volScalarField> trho
      (
        new volScalarField
        (
            IOobject
            (
                "rhoFunctionThermo",
                mesh.time().timeName(),
                this->T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimDensity
        )
    );

   //CL: copy "old" rho value onto the new rho field as start point
   //CL: for the newton solver used in this->TH( ... )
   trho()=rho_;

   volScalarField& rho = trho();

    const scalarField& hCells = h_.internalField();
    const scalarField& pCells = this->p_.internalField();
    scalarField TCells = this->T_.internalField();

    forAll(pCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        // getting the new rho Field
        mixture_.TH(hCells[celli], TCells[celli], pCells[celli], rho[celli]);
    }

    forAll(p_.boundaryField(), patchi)
    {
        fvPatchScalarField pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField ph = h_.boundaryField()[patchi];
        fvPatchScalarField pT = this->T_.boundaryField()[patchi];

        fvPatchScalarField& prho_ = rho.boundaryField()[patchi];

        forAll(pp, facei)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->patchFaceMixture(patchi, facei);

            // getting the new rho patch Field
            mixture_.TH(ph[facei], pT[facei],pp[facei],prho_[facei]);
        }
    }
    return trho;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::realGasHThermo<MixtureType>::Cv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                this->T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cv = tCv();

    forAll(this->T_, celli)
    {
        cv[celli] = this->cellMixture(celli).Cv(this->rho_[celli], this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        cv.boundaryField()[patchi] =
            Cv(this->T_.boundaryField()[patchi], patchi);
    }

    return tCv;
}


template<class MixtureType>
bool Foam::realGasHThermo<MixtureType>::read()
{
    if (basicPsiThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
