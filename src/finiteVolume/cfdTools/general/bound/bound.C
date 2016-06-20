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

#include "bound.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::bound(volScalarField& vsf, const dimensionedScalar& vsf0)
{
    scalar minVsf = min(vsf).value();

    if (minVsf < vsf0.value())
    {
        Info<< "bounding " << vsf.name()
            << ", min: " << minVsf
            << " max: " << max(vsf).value()
            << " average: " << gAverage(vsf.internalField())
            << endl;

        vsf.internalField() = max
        (
            max
            (
                vsf.internalField(),
                fvc::average(max(vsf, vsf0))().internalField()
                // Bug fix: was assuming bound on zero.  HJ, 25/Nov/2008
                *pos(vsf0.value() - vsf.internalField())
            ),
            vsf0.value()
        );

        vsf.correctBoundaryConditions();
        vsf.boundaryField() = max(vsf.boundaryField(), vsf0.value());
    }
}


void Foam::boundMinMax
(
    volScalarField& vsf,
    const dimensionedScalar& vsf0,
    const dimensionedScalar& vsf1
)
{
    scalar minVsf = min(vsf).value();
    scalar maxVsf = max(vsf).value();

    if (minVsf < vsf0.value() || maxVsf > vsf1.value())
    {
        Info<< "bounding " << vsf.name()
            << ", min: " << minVsf
            << " max: " << maxVsf
            << " average: " << gAverage(vsf.internalField())
            << endl;
    }

    if (minVsf < vsf0.value())
    {
        vsf.internalField() = max
        (
            max
            (
                vsf.internalField(),
                fvc::average(max(vsf, vsf0))().internalField()
                *pos(vsf0.value() - vsf.internalField())
            ),
            vsf0.value()
        );
        Info<< "new min: " << gMin(vsf.internalField()) << endl;
        vsf.correctBoundaryConditions();
        vsf.boundaryField() = max(vsf.boundaryField(), vsf0.value());
    }

    if (maxVsf > vsf1.value())
    {
        vsf.internalField() = min
        (
            min
            (
                vsf.internalField(),
                fvc::average(min(vsf, vsf1))().internalField()
                *neg(vsf1.value() - vsf.internalField())
                // This is needed when all values are above max
                // HJ, 18/Apr/2009
              + pos(vsf1.value() - vsf.internalField())*vsf1.value()
            ),
            vsf1.value()
        );
        Info<< "new max: " << gMax(vsf.internalField()) << endl;
        vsf.correctBoundaryConditions();
        vsf.boundaryField() = min(vsf.boundaryField(), vsf1.value());
    }
}


// ************************************************************************* //
