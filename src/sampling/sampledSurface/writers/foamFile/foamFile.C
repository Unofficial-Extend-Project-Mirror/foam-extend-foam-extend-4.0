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

#include "foamFile.H"
#include "fileName.H"
#include "OFstream.H"
#include "faceList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::foamFile<Type>::foamFile()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::foamFile<Type>::~foamFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::foamFile<Type>::write
(
    const fileName& samplePath,
    const fileName& timeDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const fileName& fieldName,
    const Field<Type>& values,
    const bool verbose
) const
{
    fileName planeFName(samplePath/timeDir/surfaceName);

    if (!exists(planeFName))
    {
        mkDir(planeFName);
    }

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << planeFName << endl;
    }

    // Points
    OFstream pointsFile(planeFName/"points");

    pointsFile << points;

    // Faces
    OFstream facesFile(planeFName/"faces");

    facesFile << faces;

    // Values to separate directory (e.g. "scalarField/p")

    fileName foamName(pTraits<Type>::typeName);

    fileName valuesDir(planeFName  / (foamName + Field<Type>::typeName));

    if (!exists(valuesDir))
    {
        mkDir(valuesDir);
    }

    OFstream valuesFile(valuesDir/fieldName);

    valuesFile << values;
}


// ************************************************************************* //
