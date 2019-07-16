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

\*---------------------------------------------------------------------------*/

#include "porousZone.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::porousZone::modifyDdt(fvMatrix<Type>& m) const
{
    if (porosity_ < 1)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneID_];

        forAll(cells, i)
        {
            m.diag()[cells[i]]   *= porosity_;
            m.source()[cells[i]] *= porosity_;
        }
    }
}


template<class RhoFieldType>
void Foam::porousZone::addPowerLawResistance
(
    scalarField& Udiag,
    const labelList& cells,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll (cells, i)
    {
        Udiag[cells[i]] +=
            V[cells[i]]*rho[cells[i]]*C0*pow(magSqr(U[cells[i]]), C1m1b2);
    }
}


template<class RhoFieldType>
void Foam::porousZone::addViscousInertialResistance
(
    scalarField& Udiag,
    vectorField& Usource,
    const labelList& cells,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    const tensor& D = D_.value();
    const tensor& F = F_.value();

    forAll (cells, i)
    {
        tensor dragCoeff = mu[cells[i]]*D + (rho[cells[i]]*mag(U[cells[i]]))*F;
        scalar isoDragCoeff = tr(dragCoeff);

        Udiag[cells[i]] += V[cells[i]]*isoDragCoeff;
        Usource[cells[i]] -=
            V[cells[i]]*((dragCoeff - I*isoDragCoeff) & U[cells[i]]);
    }
}


template<class RhoFieldType>
void Foam::porousZone::addPowerLawResistance
(
    tensorField& AU,
    const labelList& cells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll (cells, i)
    {
        AU[cells[i]] = AU[cells[i]]
          + I*(rho[cells[i]]*C0*pow(magSqr(U[cells[i]]), C1m1b2));
    }
}


template<class RhoFieldType>
void Foam::porousZone::addViscousInertialResistance
(
    tensorField& AU,
    const labelList& cells,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    const tensor& D = D_.value();
    const tensor& F = F_.value();

    forAll (cells, i)
    {
        AU[cells[i]] += mu[cells[i]]*D + (rho[cells[i]]*mag(U[cells[i]]))*F;
    }
}


template<class RhoFieldType>
void Foam::porousZone::addHeatSource
(
    const scalarField& Macro,
    const scalarField& posFlux,
    scalarField& QSource,
    scalarField& Qauxi,
    const scalarField& T,
    scalarField& Tauxi,
    const vectorField& U,
    const RhoFieldType& rho
) const
{
    // - Dual stream model (Single Effectiveness Model)
    const scalarField Tauxi_old = Tauxi; // create storePrevIter()
    const scalarField dT = Tauxi - T;


    const label nMacro(nCellsAuxInlet_*nVerticalCells_);
    scalarList Tmacro(nMacro, 0.0);
    scalarList dTaux(nMacro, 0.0);
    scalar QSum(0.0);

    const scalar Taux_relax(Tauxrelax_);
    const scalar c_aux(caux_);
    const scalar c_pri(1009.0);
    const scalar qm_aux(Maux_);
    const scalar qm_auxi = qm_aux/nCellsAuxInlet_;
    const scalar T_aux(Taux_);
    const scalar rho_pri(1.1021);
    const scalar Qeps(Qepsilon_);

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll(cells, i)
    {
        scalar Qcell = Qeps*rho_pri*c_pri*posFlux[cells[i]]*dT[cells[i]];

        Qauxi[cells[i]] = Qcell;                            // heat in each macro(cell)
        const int macro = Macro[cells[i]];
        dTaux[macro] = Qcell/(c_aux*qm_auxi);             // deltaTaux in each macro(cell)
        QSource[cells[i]] += Qcell/(rho_pri*c_pri);       // adding Heat to equation
        QSum += Qcell;                                    // summing for total heat of HX
    }

    reduce(dTaux, sumOp<scalarList>());
    Tmacro[0] = T_aux;
    forAll (Tmacro, i)
    {
        if (i > 0) Tmacro[i] = Tmacro[i-1] - dTaux[i-1];
        if ((i > 0) && (i % nVerticalCells_ == 0)) Tmacro[i] = T_aux;
    }

    reduce(QSum, sumOp<scalar>());
    Info << "Heat exchanger: " << name_ << endl;
    Info << "Q = " << QSum << endl;
    Info << "deltaT = " << QSum/(qm_aux*c_aux) << endl;

    forAll(cells, i)
    {
        const int macro = Macro[cells[i]];
        Tauxi[cells[i]] = Tauxi_old[cells[i]]
            // upwind scheme for Aux fluid (Tmacro = inlet temp)
          + Taux_relax*(Tmacro[macro] - Tauxi_old[cells[i]]);
    }
}


// ************************************************************************* //
