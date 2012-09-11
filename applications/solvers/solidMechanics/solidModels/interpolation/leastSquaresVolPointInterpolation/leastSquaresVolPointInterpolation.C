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

#include "leastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(leastSquaresVolPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

  void leastSquaresVolPointInterpolation::calcA(List<scalarSquareMatrix>& A) const
  {
    //Info << "leastSquaresVolPointInterpolation calcA" << endl;

    const fvMesh& mesh = mesh_;
    const pointField& points = mesh.points();
    
    //- construct 4x4 A matrix for each point
    //List<scalarSquareMatrix>& A = A_;
    
    //- populate A matrix
    forAll(points, pointi)
      {
	const labelList& pointCells = mesh.pointCells()[pointi];

	//- this component of matrix does not depend on coordinates
	A[pointi][3][3] = pointCells.size();

	//- fill the A matrices
	forAll(pointCells, pointCelli)
	  {
	    const label& celli = pointCells[pointCelli];
	    
	    const scalar& x = mesh.C()[celli].component(vector::X);
	    const scalar& y = mesh.C()[celli].component(vector::Y);
	    const scalar& z = mesh.C()[celli].component(vector::Z);
	    
	    A[pointi][0][0] += x*x;
	    A[pointi][0][1] += x*y;
	    A[pointi][0][2] += x*z;
	    A[pointi][0][3] += x;
	    
	    A[pointi][1][0] += x*y;
	    A[pointi][1][1] += y*y;
	    A[pointi][1][2] += y*z;
	    A[pointi][1][3] += y;
	    
	    A[pointi][2][0] += x*z;
	    A[pointi][2][1] += y*z;
	    A[pointi][2][2] += z*z;
	    A[pointi][2][3] += z;
	    
	    A[pointi][3][0] += x;
	    A[pointi][3][1] += y;
	    A[pointi][3][2] += z;
	    //A[pointi][3][3] = pointCells.size(); // set above
	  }
      }
    
    //- for boundary points we will include the surrounding face centres
    forAll(mesh.boundary(), patchi)
      {
	const vectorField& faceCentres = mesh.boundaryMesh()[patchi].faceCentres();
	const labelListList& pointFaces = mesh.boundaryMesh()[patchi].pointFaces();

	if(mesh.boundary()[patchi].coupled()) //- for proc boundaries
	  {
	    //- for coupled patches we will use the values at the neighbourField cell centres and we will
	    //- not use the boundary face values
	    //- neighbour cell centre are equal to the faceCell centres plus the delta vector
	    vectorField pDelta = mesh.boundary()[patchi].delta();
	    vectorField faceCellC(faceCentres.size(), vector::zero);
	    forAll(faceCentres, facei)
	      {
		label celli = mesh.boundaryMesh()[patchi].faceCells()[facei];
		faceCellC[facei] = mesh.C()[celli];
	      }
	    vectorField neiCellC = faceCellC + pDelta;

	    forAll(pointFaces, pointi)
	      {
		forAll(pointFaces[pointi], pointFacei)
		  {
		    label neiCelli = pointFaces[pointi][pointFacei];
		    const scalar& x = neiCellC[neiCelli].component(vector::X);
		    const scalar& y = neiCellC[neiCelli].component(vector::Y);
		    const scalar& z = neiCellC[neiCelli].component(vector::Z);
		    
		    label globalPointi = mesh.boundaryMesh()[patchi].meshPoints()[pointi];

		    A[globalPointi][0][0] += x*x;
		    A[globalPointi][0][1] += x*y;
		    A[globalPointi][0][2] += x*z;
		    A[globalPointi][0][3] += x;
		    
		    A[globalPointi][1][0] += x*y;
		    A[globalPointi][1][1] += y*y;
		    A[globalPointi][1][2] += y*z;
		    A[globalPointi][1][3] += y;
		    
		    A[globalPointi][2][0] += x*z;
		    A[globalPointi][2][1] += y*z;
		    A[globalPointi][2][2] += z*z;
		    A[globalPointi][2][3] += z;
		    
		    A[globalPointi][3][0] += x;
		    A[globalPointi][3][1] += y;
		    A[globalPointi][3][2] += z;
		    A[globalPointi][3][3] += 1; // = pointCells.size();
		  }
	      }
	  }
	else
	  {
	    //- each point must use at least 4 neighbouring locations otherwise A is singular
	    //- and simpleMatrix will cannot invert it
	    //- therefore empty patches values are included to make sure A is not singular 
	    forAll(pointFaces, pointi)
	      {
		label globalPointi = mesh.boundaryMesh()[patchi].meshPoints()[pointi];

		forAll(pointFaces[pointi], pointFacei)
		  {
		    //- fix: use pointFace not face philipc
		    label facei = pointFaces[pointi][pointFacei];
		    const scalar& x = faceCentres[facei].component(vector::X);
		    const scalar& y = faceCentres[facei].component(vector::Y);
		    const scalar& z = faceCentres[facei].component(vector::Z);
		    		
		    A[globalPointi][0][0] += x*x;
		    A[globalPointi][0][1] += x*y;
		    A[globalPointi][0][2] += x*z;
		    A[globalPointi][0][3] += x;
		    
		    A[globalPointi][1][0] += x*y;
		    A[globalPointi][1][1] += y*y;
		    A[globalPointi][1][2] += y*z;
		    A[globalPointi][1][3] += y;
		    
		    A[globalPointi][2][0] += x*z;
		    A[globalPointi][2][1] += y*z;
		    A[globalPointi][2][2] += z*z;
		    A[globalPointi][2][3] += z;
		    
		    A[globalPointi][3][0] += x;
		    A[globalPointi][3][1] += y;
		    A[globalPointi][3][2] += z;
		    A[globalPointi][3][3] += 1; // = pointCells.size();
		  }
	      }
	  } //- end of else
      } //- end of forAll boundary
  }


  void leastSquaresVolPointInterpolation::calcB(List<Field<vector> >& B, const GeometricField<vector, fvPatchField, volMesh>& vf) const
  {
    //Info << "leastSquaresVolPointInterpolation calcB" << endl;

    const fvMesh& mesh = mesh_;
    const pointField& points = mesh.points();
    
    //List<Field<vector> >& B = B_;

    for (direction compi = 0; compi < 3; compi++)
      {
	forAll(points, pointi)
	  {
	    const labelList& pointCells = mesh.pointCells()[pointi];
	    
	    forAll(pointCells, pointCelli)
	      {
		const label& celli = pointCells[pointCelli];
		
		const scalar& x = mesh.C()[celli].component(vector::X);
		const scalar& y = mesh.C()[celli].component(vector::Y);
		const scalar& z = mesh.C()[celli].component(vector::Z);
		
		const scalar& phiCompi = vf.internalField()[celli].component(compi);
		
		B[pointi][0].component(compi) += phiCompi*x;
		B[pointi][1].component(compi) += phiCompi*y;
		B[pointi][2].component(compi) += phiCompi*z;
		B[pointi][3].component(compi) += phiCompi;
	      }
	  }
	
	//- for boundary points we will include the surrounding face centres
	forAll(mesh.boundary(), patchi)
	  {
	    const vectorField& faceCentres = mesh.boundaryMesh()[patchi].faceCentres();
	    const labelListList& pointFaces = mesh.boundaryMesh()[patchi].pointFaces();
	    const labelList& faceCells = mesh.boundaryMesh()[patchi].faceCells();
	    
	    //- fix: do not calculate B for empty patches - philipc
	    if(mesh.boundary()[patchi].coupled())
	      {
		//- for coupled patches we will use the values at the neighbourField cell centres and we will
		//- not use the boundary face values
		//- neighbour cell centre are equal to the faceCell centres plus the delta vector
		vectorField pDelta = mesh.boundary()[patchi].delta();
		vectorField faceCellC(faceCentres.size(), vector::zero);
		forAll(faceCentres, facei)
		  {
		    label celli = mesh.boundaryMesh()[patchi].faceCells()[facei];
		    faceCellC[facei] = mesh.C()[celli];
		  }
		vectorField neiCellC = faceCellC + pDelta;
		
		vectorField phiNeiField = vf.boundaryField()[patchi].patchNeighbourField();

		forAll(pointFaces, pointi)
		  {
		    forAll(pointFaces[pointi], pointFacei)
		      {
			label neiCelli = pointFaces[pointi][pointFacei];
			const scalar& x = neiCellC[neiCelli].component(vector::X);
			const scalar& y = neiCellC[neiCelli].component(vector::Y);
			const scalar& z = neiCellC[neiCelli].component(vector::Z);
			
			label globalPointi = mesh.boundaryMesh()[patchi].meshPoints()[pointi];
			
			//- this is the value of phi at the cell centre in the neighbour (i.e. across the interface)
			scalar phiCompi = phiNeiField[neiCelli].component(compi);

			B[globalPointi][0].component(compi) += phiCompi*x;
			B[globalPointi][1].component(compi) += phiCompi*y;
			B[globalPointi][2].component(compi) += phiCompi*z;
			B[globalPointi][3].component(compi) += phiCompi;
		      }
		  }
	      }
	    else
	      {
	    //- each point must use at least 4 neighbouring locations otherwise A is singular
	    //- and simpleMatrix will cannot invert it
	    //- therefore empty patches values are included to make sure A is not singular 
		forAll(pointFaces, pointi)
		  {
		    forAll(pointFaces[pointi], pointFacei)
		      {
			//- fix: use pointFace not face philipc
			label facei = pointFaces[pointi][pointFacei];
			const scalar& x = faceCentres[facei].component(vector::X);
			const scalar& y = faceCentres[facei].component(vector::Y);
			const scalar& z = faceCentres[facei].component(vector::Z);
			
			label globalPointi = mesh.boundaryMesh()[patchi].meshPoints()[pointi];

			scalar phiCompi = 0.0;
			if(mesh.boundary()[patchi].type() == "empty")
			  {
			    //- use faceCell value for empty because empty patches do not store any values
			    const label& ci = faceCells[facei];
			    phiCompi = vf.internalField()[ci].component(compi);
			  }
			else
			  {
			    phiCompi = vf.boundaryField()[patchi][facei].component(compi);
			  }

    			B[globalPointi][0].component(compi) += phiCompi*x;
			B[globalPointi][1].component(compi) += phiCompi*y;
			B[globalPointi][2].component(compi) += phiCompi*z;
			B[globalPointi][3].component(compi) += phiCompi;
		      }
		  }
	      }
	  } //- end of forAll boundary
      } //- end of for all components 
  }
  

  void leastSquaresVolPointInterpolation::interpolate
  (
   const GeometricField<vector, fvPatchField, volMesh>& vf,
   GeometricField<vector, pointPatchField, pointMesh>& pf //Field<vector>& pf
   ) const
  {
    //Info << "Interpolating cell to point using leastSquaresVolPointInterpolation" << endl;

    const fvMesh& mesh = mesh_;
    const pointField& points = mesh.points();

    //- first check that point field is the correct size
    if(pf.size() != points.size())
      {
	FatalError << "pointfield should be equal to the number of points in the fvMesh" << endl;
      }

    //- calculate A and B vector
    // List<Field<vector> >& B = B_;
    // calcB(vf);
    // const List<scalarSquareMatrix>& A = A_;
    List<scalarSquareMatrix> A(mesh.points().size(), scalarSquareMatrix(4, 0.0));
    calcA(A);

    List<Field<vector> > B(mesh.points().size(), Field<vector>(4, pTraits<vector>::zero));
    calcB(B, vf);

    //- solve equations for each component of each point
    forAll(points, pointi)
      {
	Field<vector>& source = B[pointi];
	simpleMatrix<vector> leastSquaresMatrix(A[pointi], source);

	//- solve using Gauss elimination or LU decomposition with pivoting
	//Field<vector> leastSquaresSol = leastSquaresMatrix.solve();
	Field<vector> leastSquaresSol = leastSquaresMatrix.LUsolve();
	
	const scalar& x = mesh.points()[pointi].component(vector::X);
	const scalar& y = mesh.points()[pointi].component(vector::Y);
	const scalar& z = mesh.points()[pointi].component(vector::Z);
	
	//- calculate phi at vertex
	for (direction compi = 0; compi < 3; compi++)
	  {	  
	    const scalar& a = leastSquaresSol[0].component(compi);
	    const scalar& b = leastSquaresSol[1].component(compi);
	    const scalar& c = leastSquaresSol[2].component(compi);
	    const scalar& d = leastSquaresSol[3].component(compi);
	    
	    pf[pointi].component(compi) = a*x + b*y + c*z + d;
	  }
      }

    //- proc patches are synchronised
    pf.correctBoundaryConditions();
  }

  
// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //
  
  leastSquaresVolPointInterpolation::leastSquaresVolPointInterpolation(const fvMesh& vm)
  :
    mesh_(vm) //,
    //A_(vm.points().size(), scalarSquareMatrix(4, 0.0)),
    //B_(vm.points().size(), Field<vector>(4, vector::zero))
  {
    //calcA();
  }
    
// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

  leastSquaresVolPointInterpolation::~leastSquaresVolPointInterpolation()
  {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
