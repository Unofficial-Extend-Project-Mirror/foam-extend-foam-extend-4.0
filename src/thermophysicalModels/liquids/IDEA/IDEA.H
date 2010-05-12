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

Class
    Foam::IDEA

Description
    The IDEA fuel is constructed by adding 30% alphaMethylNaphthalene
    with 70% n-decane.

    The new properties have been calculated by adding the values in these
    proportions and making a least square fit, using the same NSRDS-eq.
    so that Y = 0.3*Y_naphthalene + 0.7*Y_decane

    The valid Temperature range for n-decane is normally 243.51 - 617.70 K
    and for the naphthalene it is 242.67 - 772.04 K
    The least square fit was done in the interval 244 - 617 K

    The critical temperature was taken to be 618.074 K, since this
    is the 'c'-value in the rho-equation, which corresponds to Tcrit,
    This value was then used in the fit for the NSRDS6-eq, which uses Tcrit.
    (important for the latent heat and surface tension)

    The molecular weights are 142.20 and 142.285 and for the IDEA fuel
    it is thus 142.26 ( approximately 0.3*142.2 + 0.7*142.285 )

    Critical pressure was set to the lowest one (n-Decane)

    Critical volume... also the lowest one (naphthalene) 0.523 m^3/kmol

    Second Virial Coefficient is n-Decane

SourceFiles
    IDEA.C

\*---------------------------------------------------------------------------*/

#ifndef IDEA_H
#define IDEA_H

#include "liquid.H"
#include "NSRDSfunc0.H"
#include "NSRDSfunc1.H"
#include "NSRDSfunc2.H"
#include "NSRDSfunc3.H"
#include "NSRDSfunc4.H"
#include "NSRDSfunc5.H"
#include "NSRDSfunc6.H"
#include "NSRDSfunc7.H"
#include "APIdiffCoefFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IDEA Declaration
\*---------------------------------------------------------------------------*/

class IDEA
:
    public liquid
{
    // Private data

        NSRDSfunc5 rho_;
        NSRDSfunc1 pv_;
        NSRDSfunc6 hl_;
        NSRDSfunc0 cp_;
        NSRDSfunc0 h_;
        NSRDSfunc7 cpg_;
        NSRDSfunc4 B_;
        NSRDSfunc1 mu_;
        NSRDSfunc2 mug_;
        NSRDSfunc0 K_;
        NSRDSfunc2 Kg_;
        NSRDSfunc6 sigma_;
        APIdiffCoefFunc D_;


public:

    //- Runtime type information
    TypeName("IDEA");


    // Constructors

        //- Construct null
        IDEA()
        :
            liquid(142.26, 618.074, 2.11e+6, 0.523, 0.247, 242.67, 3.4929e-2, 447.3, 1.7012e-30, 0.3478, 1.57e+4),
            rho_(152.012105, 3.87150382e-1, 618.073893, 4.00790044e-1),
            pv_(8.4817774623e+01, -8.6782398353e+03, -9.1277694857, 4.6153144498e-06, 2.0),
            hl_(618.074,2.1671983789e+05, -4.2413153435e+00, 1.1656811532e+01, -1.1656446689e+01, 4.3667661492),
            cp_(1.6604957e+3, -6.250871e-1, 6.1778552e-3, 0.0, 0.0, 0.0),
            // NN: enthalpy, h_, is not used in the sprayModel.
            // For consistency, the enthalpy is derived from hlat and hl.
            // For this fuel I've put it to zero because I was lazy at the time...
            h_(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            cpg_(1.0457515243e+03, 3.4410492875e+03, 1.5976862298e+03, 2.4697705752e+03, 7.3699710536e+02),
            B_(0.00337351091119935, -4.13606494008504, -534560.916470464, -1.13364022911762e+19, 2.80704220402713e+21),
            mu_(-6.9645853822e+01, 4.4390635942e+03, 8.4680722718e+00, 0.0, 0.0),
            mug_(4.2629382158e-08, 8.8144402122e-01, 9.6918097636e+01, 0.0),
            K_(2.03684e-01, -2.3168e-04, 0.0, 0.0, 0.0, 0.0),
            Kg_(-5.664925956707e+02, 8.896721676320e-01, -2.849783998688e+09, 6.914935658053e+05),
            sigma_(618.074, 8.3846525429e-03, -1.0044759047e+01, 2.7261918781e+01, -2.5529134309e+01, 8.6488806234),
            D_(147.18, 20.1, 142.2, 28) // NN: Same as nHeptane
        {}
        IDEA
        (
            const liquid& l,
            const NSRDSfunc5& density,
            const NSRDSfunc1& vapourPressure,
            const NSRDSfunc6& heatOfVapourisation,
            const NSRDSfunc0& heatCapacity,
            const NSRDSfunc0& enthalpy,
            const NSRDSfunc7& idealGasHeatCapacity,
            const NSRDSfunc4& secondVirialCoeff,
            const NSRDSfunc1& dynamicViscosity,
            const NSRDSfunc2& vapourDynamicViscosity,
            const NSRDSfunc0& thermalConductivity,
            const NSRDSfunc2& vapourThermalConductivity,
            const NSRDSfunc6& surfaceTension,
            const APIdiffCoefFunc& vapourDiffussivity
        )
        :
            liquid(l),
            rho_(density),
            pv_(vapourPressure),
            hl_(heatOfVapourisation),
            cp_(heatCapacity),
            h_(enthalpy),
            cpg_(idealGasHeatCapacity),
            B_(secondVirialCoeff),
            mu_(dynamicViscosity),
            mug_(vapourDynamicViscosity),
            K_(thermalConductivity),
            Kg_(vapourThermalConductivity),
            sigma_(surfaceTension),
            D_(vapourDiffussivity)
        {}

        //- Construct from Istream
        IDEA(Istream& is)
        :
            liquid(is),
            rho_(is),
            pv_(is),
            hl_(is),
            cp_(is),
            h_(is),
            cpg_(is),
            B_(is),
            mu_(is),
            mug_(is),
            K_(is),
            Kg_(is),
            sigma_(is),
            D_(is)
        {}


    // Member Functions

        //- Liquid density [kg/m^3]
        scalar rho(scalar p, scalar T) const
        {
            return rho_.f(p, T);
        }

        //- Vapour pressure [Pa]
        scalar pv(scalar p, scalar T) const
        {
            return pv_.f(p, T);
        }

        //- Heat of vapourisation [J/kg]
        scalar hl(scalar p, scalar T) const
        {
            return hl_.f(p, T);
        }

        //- Liquid heat capacity [J/(kg K)]
        scalar cp(scalar p, scalar T) const
        {
            return cp_.f(p, T);
        }

        //- Liquid Enthalpy [J/(kg)]
        scalar h(scalar p, scalar T) const
        {
            return h_.f(p, T);
        }

        //- Ideal gas heat capacity [J/(kg K)]
        scalar cpg(scalar p, scalar T) const
        {
            return cpg_.f(p, T);
        }

        //- Second Virial Coefficient [m^3/kg]
        scalar B(scalar p, scalar T) const
        {
            return B_.f(p, T);
        }

        //- Liquid viscosity [Pa s]
        scalar mu(scalar p, scalar T) const
        {
            return mu_.f(p, T);
        }

        //- Vapour viscosity [Pa s]
        scalar mug(scalar p, scalar T) const
        {
            return mug_.f(p, T);
        }

        //- Liquid thermal conductivity  [W/(m K)]
        scalar K(scalar p, scalar T) const
        {
            return K_.f(p, T);
        }

        //- Vapour thermal conductivity  [W/(m K)]
        scalar Kg(scalar p, scalar T) const
        {
            return Kg_.f(p, T);
        }

        //- Surface tension [N/m]
        scalar sigma(scalar p, scalar T) const
        {
            return sigma_.f(p, T);
        }

        //- Vapour diffussivity [m2/s]
        scalar D(scalar p, scalar T) const
        {
            return D_.f(p, T);
        }


        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            liquid::writeData(os); os << nl;
            rho_.writeData(os); os << nl;
            pv_.writeData(os); os << nl;
            hl_.writeData(os); os << nl;
            cp_.writeData(os); os << nl;
            cpg_.writeData(os); os << nl;
            mu_.writeData(os); os << nl;
            mug_.writeData(os); os << nl;
            K_.writeData(os); os << nl;
            Kg_.writeData(os); os << nl;
            sigma_.writeData(os); os << nl;
            D_.writeData(os); os << endl;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const IDEA& l)
        {
            l.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
