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

Info << nl << "Creating autocorrelation functions." << endl;

IOdictionary mdTransportProperitesDict
(
    IOobject
    (
        "mdTransportProperitesDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

const dictionary& autocorrelationFunctionDict
(
    mdTransportProperitesDict.subDict("autocorrelationFunctions")
);

//- Velocity autocorrelation function

Info << tab << "velocty" << endl;

const dictionary& velocityACFDict
(
    autocorrelationFunctionDict.subDict("velocity")
);

correlationFunction<vector> vacf
(
    mesh,
    velocityACFDict,
    molecules.size()
);

bool writeVacf(Switch(velocityACFDict.lookup("writeFile")));

//- Pressure autocorrelation function

Info << tab << "pressure" << endl;

const dictionary& pressureACFDict
(
    autocorrelationFunctionDict.subDict("pressure")
);

correlationFunction<vector> pacf
(
    mesh,
    pressureACFDict,
    1
);

bool writePacf(Switch(pressureACFDict.lookup("writeFile")));

//- Heat flux autocorrelation function

Info << tab << "heat flux" << endl;

const dictionary& heatFluxACFDict
(
    autocorrelationFunctionDict.subDict("heatFlux")
);

correlationFunction<vector> hfacf
(
    mesh,
    heatFluxACFDict,
    1
);

bool writeHFacf(Switch(heatFluxACFDict.lookup("writeFile")));
