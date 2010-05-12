// The FOAM Project // File: iC3H8O/iC3H8O.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   iC3H8O
   \\  /           | Family: iC3H8O
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2000 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    iC3H8O

DESCRIPTION
    iso-propanol


*/
// ------------------------------------------------------------------------- //

#ifndef iC3H8O_H
#define iC3H8O_H

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class iC3H8O Declaration
\*---------------------------------------------------------------------------*/

class iC3H8O
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
        NSRDSfunc0 sigma_;
        NSRDSfunc1 D_;

public:

    //- Runtime type information
    TypeName("iC3H8O");


    // Constructors

        //- Construct null
        iC3H8O()
        :
            liquid(60.096, 508.31, 4.7643e+6, 0.22013, 0.248, 185.28, 3.20e-2, 355.41, 5.5372e-30, 0.6689, 2.3575e+4),
            rho_(70.91328, 0.26475, 508.31, 0.243),
            pv_(92.935, -8177.1, -10.031, 3.9988e-06, 2),
            hl_(508.31, 948149.627263046, 0.087, 0.3007, 0, 0),
            cp_(7760.91586794462, -68.3672790202343, 0.241380457933972, -0.000235057241746539, 0, 0),
            // NN: enthalpy, h_, is not used in the sprayModel.
            // For consistency, the enthalpy is derived from hlat and hl.
            // It is, however, convenient to have it available.
            h_(-6227786.27583977, 7760.91586794462, -34.1836395101172, 0.0804601526446574, -5.87643104366347e-05, 0),
            cpg_(789.73642172524, 3219.8482428115, 1124, 1560.83599574015, 460),
            B_(0.000502529286474973, -0.104665867944622, -717185.83599574, 3.3047124600639e+18, -1.43270766773163e+21),
            mu_(-8.23, 2282.2, -0.98495, 0, 0),
            mug_(1.993e-07, 0.7233, 178, 0),
            K_(0.2029, -0.0002278, 0, 0, 0, 0),
            Kg_(-80.642, -1.4549, -604.42, 0),
            sigma_(0.03818, -3.818e-05, -6.51e-08, 0, 0, 0),
            D_(4.75e-10, 1.75, 0.0, 0.0, 0.0) // NN. same as iC3H8O
        {}
        iC3H8O
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
            const NSRDSfunc0& surfaceTension,
            const NSRDSfunc1& vapourDiffussivity
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
        iC3H8O(Istream& is)
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

        friend Ostream& operator<<(Ostream& os, const iC3H8O& l)
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

