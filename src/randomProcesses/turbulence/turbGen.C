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

\*---------------------------------------------------------------------------*/

#include "fft.H"
#include "turbGen.H"
#include "Kmesh.H"
#include "primitiveFields.H"
#include "Ek.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
turbGen::turbGen(const Kmesh& k, const scalar EA, const scalar K0)
:
    K(k),
    Ea(EA),
    k0(K0),
    RanGen(label(0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// generate and return a velocity field
vectorField turbGen::U()
{
    scalar pi2=2*mathematicalConstant::pi;

    vectorField s(K.size());
    scalarField rndPhases(K.size());

    forAll(K, i)
    {
        s[i] = RanGen.vector01();
        rndPhases[i] = RanGen.scalar01();
    }

    s = K ^ s;
    s = s/(mag(s) + 1.0e-20);

    s = Ek(Ea, k0, mag(K))*s;

    complexVectorField up = fft::reverseTransform
    (
        ComplexField(cos(pi2*rndPhases)*s, sin(pi2*rndPhases)*s),
        K.nn()
    );

    return ReImSum(up);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

