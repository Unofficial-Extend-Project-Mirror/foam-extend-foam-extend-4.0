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

Application
    transformPoints

Description
    Transforms the mesh points in the polyMesh directory according to the
    options:

    -translate vector
        Translates the points by the given vector,

    -rotate (vector vector)
        Rotates the points from the first vector to the second,

    -scale vector
        Scales the points by the given vector.

    The any or all of the three options may be specified and are processed
    in the above order.

    With -rotateFields (in combination with -rotate) it will also
    read & transform vector & tensor fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IStringStream.H"
#include "RodriguesRotation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void readAndRotateFields
(
    PtrList<GeoField>& flds,
    const fvMesh& mesh,
    const tensor& T,
    const IOobjectList& objects
)
{
    ReadFields(mesh, objects, flds);
    forAll(flds, i)
    {
        Info<< "Transforming " << flds[i].name() << endl;
        dimensionedTensor dimT("t", flds[i].dimensions(), T);
        transform(flds[i], dimT, flds[i]);
    }
}


void rotateFields(const Time& runTime, const tensor& T)
{
#   include "createMesh.H"

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    readAndRotateFields(vsFlds, mesh, T, objects);

    PtrList<volVectorField> vvFlds;
    readAndRotateFields(vvFlds, mesh, T, objects);

    PtrList<volSphericalTensorField> vstFlds;
    readAndRotateFields(vstFlds, mesh, T, objects);

    PtrList<volSymmTensorField> vsymtFlds;
    readAndRotateFields(vsymtFlds, mesh, T, objects);

    PtrList<volTensorField> vtFlds;
    readAndRotateFields(vtFlds, mesh, T, objects);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    readAndRotateFields(ssFlds, mesh, T, objects);

    PtrList<surfaceVectorField> svFlds;
    readAndRotateFields(svFlds, mesh, T, objects);

    PtrList<surfaceSphericalTensorField> sstFlds;
    readAndRotateFields(sstFlds, mesh, T, objects);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    readAndRotateFields(ssymtFlds, mesh, T, objects);

    PtrList<surfaceTensorField> stFlds;
    readAndRotateFields(stFlds, mesh, T, objects);

    mesh.write();
}


//  Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("translate", "vector");
    argList::validOptions.insert("rotate", "(vector vector)");
    argList::validOptions.insert("rotateAlongVector", "(vector angleInDegree)");
    argList::validOptions.insert("rotateFields", "");
    argList::validOptions.insert("scale", "vector");

#   include "setRootCase.H"
#   include "createTime.H"

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );


    if (args.options().size() == 0)
    {
        FatalErrorIn(args.executable())
            << "No options supplied, please use one or more of "
               "-translate, -rotate or -scale options."
            << exit(FatalError);
    }

    if (args.options().found("translate"))
    {
        vector transVector(IStringStream(args.options()["translate"])());

        Info<< "Translating points by " << transVector << endl;

        points += transVector;
    }

    if (args.options().found("rotate"))
    {
        Pair<vector> n1n2(IStringStream(args.options()["rotate"])());
        n1n2[0] /= mag(n1n2[0]);
        n1n2[1] /= mag(n1n2[1]);
        tensor T = rotationTensor(n1n2[0], n1n2[1]);

        Info<< "Rotating points by " << T << endl;

        points = transform(T, points);

        if (args.options().found("rotateFields"))
        {
            rotateFields(runTime, T);
        }
    }

    if (args.options().found("rotateAlongVector"))
    {
        IStringStream rotateVectorOptions(args.options()["rotateAlongVector"]);

        vector rotationAxis(rotateVectorOptions);
        scalar rotationAngle = readScalar(rotateVectorOptions);

        tensor T = RodriguesRotation(rotationAxis, rotationAngle);

        Info << "Rotating points by " << T << endl;
 
        points = transform(T, points);

        if (args.options().found("rotateFields"))
        {
            rotateFields(runTime, T);
        }
    }

    if (args.options().found("scale"))
    {
        vector scaleVector(IStringStream(args.options()["scale"])());

        Info<< "Scaling points by " << scaleVector << endl;

        points.replace(vector::X, scaleVector.x()*points.component(vector::X));
        points.replace(vector::Y, scaleVector.y()*points.component(vector::Y));
        points.replace(vector::Z, scaleVector.z()*points.component(vector::Z));
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << "Writing points into directory " << points.path() << nl << endl;
    points.write();

    return(0);
}


// ************************************************************************* //
