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
    Foam::IC8H18

Description
    iso-Octane

SourceFiles
    IC8H18.C

\*---------------------------------------------------------------------------*/

#ifndef IC8H18_H
#define IC8H18_H

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
                           Class IC8H18 Declaration
\*---------------------------------------------------------------------------*/

class IC8H18
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
    TypeName("IC8H18");


    // Constructors

        //- Construct null
        IC8H18()
        :
            liquid(114.231, 543.96, 2.5676e+6, 0.468, 0.266, 165.78, 1.4464e-2, 372.39, 0.0, 0.3031, 1.4051e+4),
            rho_(67.2363666, 0.27373, 543.96, 0.2846),
            pv_(120.81, -7550, -16.111, 0.017099, 1),
            hl_(543.96, 375379.713037617, 0.1549, 0.138, 0.0666, 0),
            cp_(1219.89652546156, 1.67205049417409, 0.00414073237562483, 0, 0, 0),
            // NN: enthalpy, h_, is not used in the sprayModel.
            // For consistency, the enthalpy is derived from hlat and hl.
            // It is, however, convenient to have it available.
            h_(-2743275.10767575, 1219.89652546156, 0.836025247087043, 0.00138024412520828, 0, 0),
            cpg_(997.10236275617, 4627.4653990598, 1594, 2933.52942721328, 677.94),
            B_(0.00234936225718063, -2.83381919093766, -413154.047500241, -3.49703670632315e+17, 3.13750207912038e+19),
            mu_(-15.811, 1282.5, 0.67791, -3.8617e-28, 10),
            mug_(1.107e-07, 0.746, 72.4, 0),
            K_(0.1508, -0.0001712, 0, 0, 0, 0),
            Kg_(1.758e-05, 1.3114, 392.9, 0),
            sigma_(543.96, 0.047434, 1.1975, 0, 0, 0),
            D_(147.18, 20.1, 114.231, 28) // NN: Same as nHeptane
        {}
        IC8H18
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
        IC8H18(Istream& is)
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

        friend Ostream& operator<<(Ostream& os, const IC8H18& l)
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
