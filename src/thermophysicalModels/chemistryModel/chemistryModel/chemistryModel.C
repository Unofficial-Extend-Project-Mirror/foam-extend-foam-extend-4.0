/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "chemistryModel.H"
#include "chemistrySolver.H"
#include "multiComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::chemistryModel::chemistryModel
(
    hCombustionThermo& thermo,
    const volScalarField& rho
)
:
    thermo_(thermo),

    Y_(thermo.composition().Y()),
    rho_(rho),

    chemistryProperties_
    (
        IOobject
        (
            "chemistryProperties",
            rho_.time().constant(),
            rho_.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    chemistry_(chemistryProperties_.lookup("chemistry")),

    reactions_(dynamic_cast<const reactingMixture&>(thermo)),
    specieThermo_(dynamic_cast<const reactingMixture&>(thermo).speciesData()),

    Ns_(thermo.composition().Y().size()),
    Nr_(reactions_.size()),
    coeffs_(Ns_ + 2, 0.0),

    solver_
    (
        chemistrySolver::New
        (
            chemistryProperties_,
            *this
        )
    ),
    deltaTChem_
    (
        rho_.size(),
        readScalar(chemistryProperties_.lookup("initialChemicalTimeStep"))
    )
{
    // set the size of the chemistry sources
    RR_.setSize(Ns());
    for(label i=0; i<Ns(); i++)
    {
        RR_.set
        (
            i,
            new scalarField(rho_.size(), 0.0)
        );
    }

    Info<< "chemistryModel::chemistryModel: Number of species = " << Ns()
        << " and reactions = " << Nr() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::chemistryModel::~chemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::chemistryModel::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    scalar pf,cf,pr,cr;
    label lRef, rRef;

    scalarField om(nEqns(), 0.0);

    forAll(reactions_, i)
    {
        const reaction& R = reactions_[i];
        
        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            label si = R.lhs()[s].index;
            scalar sl = R.lhs()[s].stoichCoeff;
            om[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            label si = R.rhs()[s].index;
            scalar sr = R.rhs()[s].stoichCoeff;
            om[si] += sr*omegai;
        }
    }

    return om;
} // end omega


Foam::scalar Foam::chemistryModel::omega
(
    const reaction& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalarField c2(Ns(), 0.0);
    for(label i=0; i<Ns(); i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    scalar kf = R.kf(T, p, c2);
    scalar kr = R.kr(kf, T, p, c2);

    pf = 1.0;
    pr = 1.0;
    
    label Nl = R.lhs().size();
    label Nr = R.rhs().size();
    
    label slRef = 0;
    lRef = R.lhs()[slRef].index;
        
    pf = kf;
    for(label s=1; s<Nl; s++)
    {
        label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        scalar exp = R.lhs()[slRef].exponent;
        if (exp<1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp-1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp-1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;
    
    // find the matrix element and element position for the rhs
    pr = kr;
    for(label s=1; s<Nr; s++)
    {
        label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        scalar exp = R.rhs()[srRef].exponent;
        if (exp<1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp-1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp-1.0);
        }
        
    }

    return pf*cf - pr*cr;

} // end omega


void Foam::chemistryModel::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{
    scalar T = c[Ns_];
    scalar p = c[Ns_ + 1];

    dcdt = omega(c, T, p);

    // constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
    scalar cSum = 0.0;
    for(label i=0; i<Ns(); i++)
    {
        scalar W = specieThermo_[i].W();
        cSum += c[i];
        rho += W*c[i];
    }
    scalar mw = rho/cSum;
    scalar cp = 0.0;
    for(label i=0; i<Ns(); i++)
    {
        scalar cpi = specieThermo_[i].cp(T);
        scalar Xi = c[i]/rho;
        cp += Xi*cpi;
    }
    cp /= mw;

    scalar dT = 0.0;
    for(label i=0; i<Ns(); i++)
    {
        scalar hi = specieThermo_[i].h(T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    // limit the time-derivative, this is more stable for the ODE
    // solver when calculating the allowed time step
    scalar dtMag = min(500.0, mag(dT));
    dcdt[Ns_] = -dT*dtMag/(mag(dT) + 1.0e-10);

    // dp/dt = ...
    dcdt[Ns_+1] = 0.0;
}


void Foam::chemistryModel::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    Matrix<scalar>& dfdc
) const
{
    scalar T = c[Ns_];
    scalar p = c[Ns_ + 1];
    
    scalarField c2(Ns(), 0.0);
    for(label i=0; i<Ns(); i++)
    {
        c2[i] = max(c[i], 0.0);
    }

    for(label i=0; i<nEqns(); i++)
    {
        for(label j=0; j<nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    // length of the first argument must be Ns()
    dcdt = omega(c2, T, p);

    for (label ri=0; ri<reactions_.size(); ri++)
    {
        const reaction& R = reactions_[ri];

        scalar kf0 = R.kf(T, p, c2);
        scalar kr0 = R.kr(T, p, c2);

        forAll(R.lhs(), j)
        {
            label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar el = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (el < 1.0)
                    {
                        if (c2[si]>SMALL)
                        {
                            kf *= el*pow(c2[si]+VSMALL, el-1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*pow(c2[si], el-1.0);
                    }
                }
                else
                {
                    kf *= pow(c2[si], el);
                }
            }

            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] -= sl*kf;
            }
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] += sr*kf;
            }
        }

        forAll(R.rhs(), j)
        {
            label sj = R.rhs()[j].index;
            scalar kr = kr0;
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar er = R.rhs()[i].exponent;
                if (i==j)
                {
                    if (er<1.0)
                    {
                        if (c2[si]>SMALL)
                        {
                            kr *= er*pow(c2[si]+VSMALL, er-1.0);
                        }
                        else
                        {
                            kr = 0.0;
                        }
                    }
                    else
                    {
                        kr *= er*pow(c2[si], er-1.0);
                    }
                }
                else
                {
                    kr *= pow(c2[si], er);
                }
            }

            forAll(R.lhs(), i)
            {
                label si = R.lhs()[i].index;
                scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] += sl*kr;
            }
            forAll(R.rhs(), i)
            {
                label si = R.rhs()[i].index;
                scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] -= sr*kr;
            }
        }
    }

    // calculate the dcdT elements numerically
    scalar delta = 1.0e-8;
    scalarField dcdT0 = omega(c2, T-delta, p);
    scalarField dcdT1 = omega(c2, T+delta, p);

    for(label i=0; i<nEqns(); i++)
    {
        dfdc[i][Ns()] = 0.5*(dcdT1[i]-dcdT0[i])/delta;
    }

} // end jacobian


Foam::tmp<Foam::volScalarField> Foam::chemistryModel::tc() const
{

    scalar pf,cf,pr,cr;
    label lRef, rRef;

    label nCells = rho_.size();
    label Nr = reactions_.size();

    scalarField t(nCells, SMALL);

    if (chemistry_)
    {

        for(label celli=0; celli<nCells; celli++)
        {
            scalar rhoi = rho_[celli];
            scalar Ti = thermo_.T()[celli];
            scalar pi = thermo_.p()[celli];
            scalarField c(Ns_);
            scalar cSum = 0.0;
            
            for(label i=0; i<Ns_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
                cSum += c[i];
            }
            
            forAll(reactions_, i)
            {
                const reaction& R = reactions_[i];
                
                omega
                (
                    R, c, Ti, pi, pf, cf, lRef, pr, cr, rRef
                );
                
                forAll(R.rhs(), s)
                {
                    scalar sr = R.rhs()[s].stoichCoeff;
                    t[celli] += sr*pf*cf;
                }
            }
            t[celli] = Nr*cSum/t[celli];
        }
    }

    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                rho_.time().timeName(),
                rho_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho_.mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tsource().internalField() = t;
    tsource().correctBoundaryConditions();

    return tsource;
}


Foam::label Foam::chemistryModel::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return Ns_ + 2;
}


void Foam::chemistryModel::calculate()
{
    for(label i=0; i<Ns(); i++)
    {
        RR_[i].setSize(rho_.size());
    }

    if (chemistry_)
    {
        for(label celli=0; celli<rho_.size(); celli++)
        {
            for(label i=0; i<Ns(); i++)
            {
                RR_[i][celli] = 0.0;
            }
            
            scalar rhoi = rho_[celli];
            scalar Ti = thermo_.T()[celli];
            scalar pi = thermo_.p()[celli];
            
            scalarField c(Ns_);
            scalarField dcdt(nEqns(), 0.0);
            
            for(label i=0; i<Ns_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
            }
            
            dcdt = omega(c, Ti, pi);
            
            for(label i=0; i<Ns_; i++)
            {
                RR_[i][celli] = dcdt[i]*specieThermo_[i].W();
            }
        }
    }
}


Foam::scalar Foam::chemistryModel::solve(const scalar t0, const scalar deltaT)
{
    for(label i=0; i<Ns(); i++)
    {
        RR_[i].setSize(rho_.size());
    }

    if (!chemistry_)
    {
        return GREAT;
    }

    scalar deltaTMin = GREAT;

    for(label celli=0; celli<rho_.size(); celli++)
    {
        for(label i=0; i<Ns(); i++)
        {
            RR_[i][celli] = 0.0;
        }

        scalar rhoi = rho_[celli];
        scalar Ti = thermo_.T()[celli];
        scalar hi = thermo_.h()[celli];
        scalar pi = thermo_.p()[celli];

        scalarField c(Ns_);
        scalarField c0(Ns_);
        scalarField dc(Ns_, 0.0);

        for(label i=0; i<Ns_; i++)
        {
            c[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
        }
        c0 = c;

        scalar t = t0;
        scalar tauC = deltaTChem_[celli];
        scalar dt = min(deltaT, tauC);
        scalar timeLeft = deltaT;

        // calculate the chemical source terms
        scalar cTot = 0.0;

        while(timeLeft > SMALL)
        {
            tauC = solver().solve(c, Ti, pi, t, dt);
            t += dt;

            // update the temperature
            cTot = sum(c);
            reactionThermo mixture(0.0*specieThermo_[0]);
            for(label i=0; i<Ns_; i++)
            {
                mixture += (c[i]/cTot)*specieThermo_[i];
            }        
            Ti = mixture.TH(hi, Ti);

            timeLeft -= dt;
            deltaTChem_[celli] = tauC;
            dt = min(timeLeft, tauC);
            dt = max(dt, SMALL);
        }
        deltaTMin = min(tauC, deltaTMin);

        dc = c - c0;
        scalar WTot = 0.0;
        for(label i=0; i<Ns_; i++)
        {
            WTot += c[i]*specieThermo_[i].W();
        }        
        WTot /= cTot;

        for(label i=0; i<Ns_; i++)
        {
            RR_[i][celli] = dc[i]*specieThermo_[i].W()/deltaT;
        }
    }

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*deltaT);

    return deltaTMin;
}


// ************************************************************************* //
