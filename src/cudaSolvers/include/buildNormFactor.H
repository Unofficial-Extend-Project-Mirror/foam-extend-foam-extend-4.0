/**********************************************************************\
  ______  __    __   _______  _______  __       __   __   __   __  ___
 /      ||  |  |  | |   ____||   ____||  |     |  | |  \ |  | |  |/  /
|  ,----'|  |  |  | |  |__   |  |__   |  |     |  | |   \|  | |  '  /
|  |     |  |  |  | |   __|  |   __|  |  |     |  | |  . `  | |    <
|  `----.|  `--'  | |  |     |  |     |  `----.|  | |  |\   | |  .  \
 \______| \______/  |__|     |__|     |_______||__| |__| \__| |__|\__\

Cuda For FOAM Link

cufflink is a library for linking numerical methods based on Nvidia's
Compute Unified Device Architecture (CUDA™) C/C++ programming language
and OpenFOAM®.

Please note that cufflink is not approved or endorsed by ESI-OpenCFD®
Limited, the owner of the OpenFOAM® and OpenCFD® trademarks and
producer of OpenFOAM® software.

The official web-site of OpenCFD® Limited is www.openfoam.com .

------------------------------------------------------------------------
This file is part of cufflink.

    cufflink is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cufflink is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cufflink.  If not, see <http://www.gnu.org/licenses/>.

    Author
    Daniel P. Combest.  All rights reserved.
    Modifications by Dominik Christ, Wikki Ltd.

    Description
    diagonal preconditioned conjugate gradient
    solver for symmetric Matrices using a CUSP CUDA™ based solver.

\**********************************************************************/


// Build the normfactor as foam-extend does in file
// ($FOAM_SRC/foam/matrices/lduMatrix/lduMatrix/lduMatrixSolver.C)

//all variables will be deleted once out of scope
{

    //compute average of x vector
    ValueType xRef = thrust::reduce
    (
        X.begin(),
        X.begin() + A.num_rows
    )/ValueType(A.num_rows);

    //vector of average X values
    cusp::array1d<ValueType,MemorySpace> xBar(N,xRef);

    //holds A*xBar
    cusp::array1d<ValueType,MemorySpace> yBar(N);

    //compute ybar<-A*xBar
    cusp::multiply(A, xBar, yBar);

    cusp::array1d<ValueType,MemorySpace> ymyBar(N);  // holds y - yBar
    cusp::array1d<ValueType,MemorySpace> bmyBar(N);  // holds b - yBar

    // Calculate: ymyBar <- y - yBar
    cusp::blas::axpby(y, yBar, ymyBar, ValueType(1), ValueType(-1));
    // Calculate: bmyBar <- b - yBar
    cusp::blas::axpby(B, yBar, bmyBar, ValueType(1), ValueType(-1));

    //compute norm factor exactly as OpenFOAM does
    normFactor = cusp::blas::nrm2(ymyBar) + cusp::blas::nrm2(bmyBar) + SMALL;

    if (solverPerf->debugCusp)
    {
        std::cout<<"   Normalisation factor = "<<normFactor<<"\n";
    }
}
