/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "volPointInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "demandDrivenData.H"
#include "coupledPointPatchFields.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void volPointInterpolation::makeWeights() const
{
    if (debug)
    {
        Info<< "volPointInterpolation::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    // It is an error to attempt to recalculate if the pointer is already set
    if (pointWeightsPtr_)
    {
        FatalErrorIn("volPointInterpolation::makeWeights() const")
            << "Point weights already calculated."
            << abort(FatalError);
    }

    // Get mesh data
    const pointField& points = mesh().points();
    const labelListList& pointCells = mesh().pointCells();
    const vectorField& cellCentres = mesh().cellCentres();

    // Allocate storage for weighting factors
    pointWeightsPtr_ = new scalarListList(points.size());
    scalarListList& pointWeights = *pointWeightsPtr_;

    forAll(pointWeights, pointi)
    {
        pointWeights[pointi].setSize(pointCells[pointi].size());
    }

    // Note: volPointInterpolation is a mesh object which creates here another
    // mesh object of type pointMesh, which it directly depends on. This
    // function should be therefore called only when all mesh objects are
    // updated (on topo changes). In case the pointMesh is not updated and
    // someone asks for weights on the "updated" volPointInterpolation, this
    // will fail. VV, 5/Feb/2018.
    pointScalarField sumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh().polyMesh::instance(),
            mesh()
        ),
        pointMesh::New(mesh()),
        dimensionedScalar("zero", dimless, 0)
    );

    // Calculate inverse distances between cell centres and points
    // and store in weighting factor array
    forAll(points, pointi)
    {
        scalarList& pw = pointWeights[pointi];
        const labelList& pcp = pointCells[pointi];

        forAll(pcp, pointCelli)
        {
            pw[pointCelli] =
                1/mag(points[pointi] - cellCentres[pcp[pointCelli]]);

            sumWeights[pointi] += pw[pointCelli];
        }
    }

    // Coupled boundary update
    forAll(sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            sumWeights.boundaryField()[patchi].initAddField();
        }
    }

    forAll(sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            sumWeights.boundaryField()[patchi].addField
            (
                sumWeights.internalField()
            );
        }
    }

    forAll(points, pointi)
    {
        scalarList& pw = pointWeights[pointi];

        forAll(pw, pointCelli)
        {
            pw[pointCelli] /= sumWeights[pointi];
        }
    }

    if (debug)
    {
        Info<< "volPointInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

volPointInterpolation::volPointInterpolation(const fvMesh& vm)
:
    MeshObject<fvMesh, volPointInterpolation>(vm),
    boundaryInterpolator_(vm),
    pointWeightsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

volPointInterpolation::~volPointInterpolation()
{
    deleteDemandDrivenData(pointWeightsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const scalarListList& volPointInterpolation::pointWeights() const
{
    if (!pointWeightsPtr_)
    {
        makeWeights();
    }

    return *pointWeightsPtr_;
}


bool volPointInterpolation::movePoints() const
{
    deleteDemandDrivenData(pointWeightsPtr_);

    // Updated for MeshObject handling
    // HJ, 30/Aug/2010
    const_cast<pointPatchInterpolation&>(boundaryInterpolator_).movePoints();

    return true;
}


bool volPointInterpolation::updateMesh(const mapPolyMesh&) const
{
    deleteDemandDrivenData(pointWeightsPtr_);

    // Updated for MeshObject handling
    // HJ, 30/Aug/2010
    const_cast<pointPatchInterpolation&>(boundaryInterpolator_).updateMesh();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
