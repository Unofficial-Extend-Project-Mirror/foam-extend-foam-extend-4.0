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
    Foam::C7H16

Description
    nHeptane

SourceFiles
    C7H16.C

\*---------------------------------------------------------------------------*/

#ifndef C7H16_H
#define C7H16_H

#include "liquid.H"
#include "NSRDSfunc0.H"
#include "NSRDSfunc1.H"
#include "NSRDSfunc2.H"
#include "NSRDSfunc3.H"
#include "NSRDSfunc4.H"
#include "NSRDSfunc5.H"
#include "NSRDSfunc6.H"
#include "NSRDSfunc7.H"
#include "NSRDSfunc14.H"
#include "APIdiffCoefFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class C7H16 Declaration
\*---------------------------------------------------------------------------*/

class C7H16
:
    public liquid
{
    // Private data

        NSRDSfunc5 rho_;
        NSRDSfunc1 pv_;
        NSRDSfunc6 hl_;
        NSRDSfunc14 cp_;
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
    TypeName("C7H16");


    // Constructors

        //- Construct null
        C7H16()
        :
            liquid(100.204, 540.20, 2.74e+6, 0.428, 0.261, 182.57, 1.8269e-1, 371.58, 0, 0.3495, 1.52e+4),
            rho_(61.38396836, 0.26211, 540.2, 0.28141),
            pv_(87.829, -6996.4, -9.8802, 7.2099e-06, 2),
            hl_(540.20, 499121.791545248, 0.38795, 0, 0, 0),
            cp_(540.20, 6.11976102401216, 3137.69909384855, 182.274175063868, -254.530511150515),
            // NN: enthalpy, h_, is not used in the sprayModel.
            // For consistency, the enthalpy is derived from hlat and hl.
            // It is, however, convenient to have it available.
            h_(-3.1469964e+6,7.3072e+3,-3.52884e+1,1.10637e-1,-1.634831e-4,9.64941e-8),
            cpg_(1199.05392998284, 3992.85457666361, 1676.6, 2734.42177956968, 756.4),
            B_(0.00274040956448844, -2.90407568560137, -440900.562851782, -8.78208454752305e+17, 1.28238393676899e+20),
            mu_(-24.451, 1533.1, 2.0087, 0, 0),
            mug_(6.672e-08, 0.82837, 85.752, 0),
            K_(0.215, -0.000303, 0, 0, 0, 0),
            Kg_(-0.070028, 0.38068, -7049.9, -2400500),
            sigma_(540.20, 0.054143, 1.2512, 0, 0, 0),
            D_(147.18, 20.1, 100.204, 28) // NN: Same as C7H16
        {}
        C7H16
        (
            const liquid& l,
            const NSRDSfunc5& density,
            const NSRDSfunc1& vapourPressure,
            const NSRDSfunc6& heatOfVapourisation,
            const NSRDSfunc14& heatCapacity,
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
        C7H16(Istream& is)
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
            B_.writeData(os); os << nl;
            mu_.writeData(os); os << nl;
            mug_.writeData(os); os << nl;
            K_.writeData(os); os << nl;
            Kg_.writeData(os); os << nl;
            sigma_.writeData(os); os << nl;
            D_.writeData(os); os << endl;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const C7H16& l)
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
