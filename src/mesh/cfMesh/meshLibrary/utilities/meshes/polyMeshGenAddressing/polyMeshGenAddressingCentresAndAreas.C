/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description
    Calulate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcFaceCentresAndAreas() const
{
    if( faceCentresPtr_ || faceAreasPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcFaceCentresAndAreas() const")
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    faceCentresPtr_ = new vectorField(faces.size());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(faces.size());
    vectorField& fAreas = *faceAreasPtr_;

    makeFaceCentresAndAreas(points, fCtrs, fAreas);
}

void polyMeshGenAddressing::makeFaceCentresAndAreas
(
    const pointFieldPMG& p,
    vectorField& fCtrs,
    vectorField& fAreas
) const
{
    const faceListPMG& fs = mesh_.faces();
    const label nFaces = fs.size();

    # ifdef USE_OMP
    # pragma omp parallel for if( nFaces > 1000 )
    # endif
    for(label faceI=0;faceI<nFaces;++faceI)
    {
        const face& f = fs[faceI];
        label nPoints = f.size();

        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            fCtrs[faceI] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
            fAreas[faceI] = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
        }
        else
        {
            vector sumN = vector::zero;
            scalar sumA = 0.0;
            vector sumAc = vector::zero;

            point fCentre = p[f[0]];
            for(label pi=1;pi<nPoints;++pi)
            {
                fCentre += p[f[pi]];
            }

            fCentre /= nPoints;

            for(label pi=0;pi<nPoints;++pi)
            {
                const point& nextPoint = p[f.nextLabel(pi)];

                vector c = p[f[pi]] + nextPoint + fCentre;
                vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
                scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            fCtrs[faceI] = (1.0/3.0)*sumAc/(sumA + VSMALL);
            fAreas[faceI] = 0.5*sumN;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& polyMeshGenAddressing::faceCentres() const
{
    if( !faceCentresPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& polyMeshGenAddressing::faceCentres() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceCentresAndAreas();
    }

    return *faceCentresPtr_;
}

const vectorField& polyMeshGenAddressing::faceAreas() const
{
    if( !faceAreasPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& polyMeshGenAddressing::faceAreas() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceCentresAndAreas();
    }

    return *faceAreasPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
