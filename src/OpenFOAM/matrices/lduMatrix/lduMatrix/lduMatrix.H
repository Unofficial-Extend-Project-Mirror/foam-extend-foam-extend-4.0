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

Class
    Foam::lduMatrix

Description
    lduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.

    Addressing arrays must be supplied for the upper and lower triangles.

    It might be better if this class were organised as a hierachy starting
    from an empty matrix, then deriving diagonal, symmetric and asymmetric
    matrices.

SourceFiles
    lduMatrixATmul.C
    lduMatrix.C
    lduMatrixTemplates.C
    lduMatrixOperations.C
    lduMatrixSolver.C
    lduMatrixPreconditioner.C
    lduMatrixTests.C
    lduMatrixUpdateMatrixInterfaces.C
    lduMatrixBufferedUpdateMatrixInterfaces.C

\*---------------------------------------------------------------------------*/

#ifndef lduMatrix_H
#define lduMatrix_H

#include "lduMesh.H"
#include "primitiveFieldsFwd.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class lduMatrix;
Ostream& operator<<(Ostream&, const lduMatrix&);


/*---------------------------------------------------------------------------*\
                           Class lduMatrix Declaration
\*---------------------------------------------------------------------------*/

class lduMatrix
{
    // Private data

        //- LDU mesh reference
        const lduMesh& lduMesh_;

        //- Coefficients (not including interfaces)
        scalarField *lowerPtr_, *diagPtr_, *upperPtr_;


public:

    //- Class returned by the solver, containing performance statistics
    class solverPerformance
    {
        // Private Data

            //- Name of solver
            word solverName_;

            //- Name of field
            word fieldName_;

            //- Initial residual
            scalar initialResidual_;

            //- Final residual
            scalar finalResidual_;

            //- Number of iterations
            label noIterations_;

            //- Has the solver converged?
            bool converged_;

            //- Has the solver indicated a singular matrix?
            bool singular_;


    public:

        // Constructors

            //- Construct null
            solverPerformance()
            :
                initialResidual_(0),
                finalResidual_(0),
                noIterations_(0),
                converged_(false),
                singular_(false)
            {}


            //- Construct from components
            solverPerformance
            (
                const word& solverName,
                const word& fieldName,
                const scalar iRes = 0,
                const scalar fRes = 0,
                const label nIter = 0,
                const bool converged = false,
                const bool singular = false
            )
            :
                solverName_(solverName),
                fieldName_(fieldName),
                initialResidual_(iRes),
                finalResidual_(fRes),
                noIterations_(nIter),
                converged_(converged),
                singular_(singular)
            {}


        // Member functions

            //- Return solver name
            const word& solverName() const
            {
                return solverName_;
            }

            //- Return initial residual
            scalar initialResidual() const
            {
                return initialResidual_;
            }

            //- Return initial residual
            scalar& initialResidual()
            {
                return initialResidual_;
            }


            //- Return final residual
            scalar finalResidual() const
            {
                return finalResidual_;
            }

            //- Return final residual
            scalar& finalResidual()
            {
                return finalResidual_;
            }


            //- Return number of iterations
            label nIterations() const
            {
                return noIterations_;
            }

            //- Return number of iterations
            label& nIterations()
            {
                return noIterations_;
            }


            //- Has the solver converged?
            bool converged() const
            {
                return converged_;
            }

            //- Is the matrix singular?
            bool singular() const
            {
                return singular_;
            }

            //- Convergence test
            bool checkConvergence
            (
                const scalar tolerance,
                const scalar relTolerance
            );

            //- Singularity test
            bool checkSingularity(const scalar residual);

            //- Print summary of solver performance
            void print() const;
    };


    //- Abstract base-class for lduMatrix solvers
    class solver
    {
        // Private data

            //- Name of field being solved for
            word fieldName_;

            //- Control data dictionary
            dictionary dict_;

            //- Solver tolerance
            scalar tolerance_;

            //- Relative tolerance
            scalar relTolerance_;

            //- Minimum number of iterations
            //  (forced irrespective of convergence)
            label minIter_;

            //- Maximum number of iterations
            label maxIter_;


    protected:

        // Protected data

            //- Matrix reference
            const lduMatrix& matrix_;

            //- Coupling boundary coefficients
            const FieldField<Field, scalar>& coupleBouCoeffs_;

            //- Coupling internal coefficients
            const FieldField<Field, scalar>& coupleIntCoeffs_;

            //- Coupling interfaces. Holding a reference means the caller
            //  is responsible for keeping a copy
            const lduInterfaceFieldPtrsList& interfaces_;


        // Protected Member Functions

            //- Return dictionary
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Is the stop criterion reached
            bool stop(lduMatrix::solverPerformance& solverPerf) const;

            //- Read the control parameters from the dictionary
            virtual void readControls();

            //- Return the matrix norm used to normalise the residual for the
            //  stopping criterion
            scalar normFactor
            (
                const scalarField& x,
                const scalarField& b,
                const scalarField& Ax,
                scalarField& tmpField,
                const direction cmpt
            ) const;

            scalar normFactor
            (
                const scalarField& x,
                const scalarField& b,
                const direction cmpt
            ) const;


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                symMatrix,
                (
                    const word& fieldName,
                    const lduMatrix& matrix,
                    const FieldField<Field, scalar>& coupleBouCoeffs,
                    const FieldField<Field, scalar>& coupleIntCoeffs,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& dict
                ),
                (
                    fieldName,
                    matrix,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaces,
                    dict
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                asymMatrix,
                (
                    const word& fieldName,
                    const lduMatrix& matrix,
                    const FieldField<Field, scalar>& coupleBouCoeffs,
                    const FieldField<Field, scalar>& coupleIntCoeffs,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& dict
                ),
                (
                    fieldName,
                    matrix,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaces,
                    dict
                )
            );


        // Constructors

            //- Construct from a dictionary
            solver
            (
                const word& fieldName,
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& dict
            );


        // Selectors

            //- Return a new solver
            static autoPtr<solver> New
            (
                const word& fieldName,
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& dict
            );


        // Destructor

            virtual ~solver()
            {}


        // Member functions

            // Access

                const word& fieldName() const
                {
                    return fieldName_;
                }

                scalar tolerance() const
                {
                    return tolerance_;
                }

                scalar relTolerance() const
                {
                    return relTolerance_;
                }

                label minIter() const
                {
                    return minIter_;
                }

                label maxIter() const
                {
                    return maxIter_;
                }

                const lduMatrix& matrix() const
                {
                    return matrix_;
                }

                 const FieldField<Field, scalar>& coupleBouCoeffs() const
                 {
                     return coupleBouCoeffs_;
                 }

                 const FieldField<Field, scalar>& coupleIntCoeffs() const
                 {
                     return coupleIntCoeffs_;
                 }

                 const lduInterfaceFieldPtrsList& interfaces() const
                 {
                     return interfaces_;
                 }


            //- Read and reset the solver parameters from the given stream
            virtual void read(const dictionary&);

            virtual solverPerformance solve
            (
                scalarField& x,
                const scalarField& b,
                const direction cmpt = 0
            ) const = 0;
    };


    //- Abstract base-class for lduMatrix smoothers
    class smoother
    {
    protected:

        // Protected data

            //- Matrix reference
            const lduMatrix& matrix_;

            //- Coupling boundary coefficients
            const FieldField<Field, scalar>& coupleBouCoeffs_;

            //- Coupling internal coefficients
            const FieldField<Field, scalar>& coupleIntCoeffs_;

            //- Coupling interfaces
            const lduInterfaceFieldPtrsList& interfaces_;


    public:

        //- Find the smoother name (directly or from a sub-dictionary)
        static word getName(const dictionary&);

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                smoother,
                symMatrix,
                (
                    const lduMatrix& matrix,
                    const FieldField<Field, scalar>& coupleBouCoeffs,
                    const FieldField<Field, scalar>& coupleIntCoeffs,
                    const lduInterfaceFieldPtrsList& interfaces
                ),
                (
                    matrix,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaces
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                smoother,
                asymMatrix,
                (
                    const lduMatrix& matrix,
                    const FieldField<Field, scalar>& coupleBouCoeffs,
                    const FieldField<Field, scalar>& coupleIntCoeffs,
                    const lduInterfaceFieldPtrsList& interfaces
                ),
                (
                    matrix,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaces
                )
            );


        // Constructors

            smoother
            (
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces
            );


        // Selectors

            //- Return a new smoother
            static autoPtr<smoother> New
            (
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& dict
            );


        // Destructor

            virtual ~smoother()
            {}


        // Member functions

            // Access

                const lduMatrix& matrix() const
                {
                    return matrix_;
                }

                 const FieldField<Field, scalar>& coupleBouCoeffs() const
                 {
                     return coupleBouCoeffs_;
                 }

                 const FieldField<Field, scalar>& coupleIntCoeffs() const
                 {
                     return coupleIntCoeffs_;
                 }

                 const lduInterfaceFieldPtrsList& interfaces() const
                 {
                     return interfaces_;
                 }


            //- Smooth the solution for a given number of sweeps
            virtual void smooth
            (
                scalarField& x,
                const scalarField& b,
                const direction cmpt,
                const label nSweeps
            ) const = 0;
    };


    //- Abstract base-class for lduMatrix preconditioners
    class preconditioner
    {
    protected:

        // Protected data

            //- Matrix reference
            const lduMatrix& matrix_;

            //- Coupling boundary coefficients
            const FieldField<Field, scalar>& coupleBouCoeffs_;

            //- Coupling internal coefficients
            const FieldField<Field, scalar>& coupleIntCoeffs_;

            //- Coupling interfaces
            const lduInterfaceFieldPtrsList& interfaces_;


    public:

        //- Find the preconditioner name (directly or from a sub-dictionary)
        static word getName(const dictionary&);

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                symMatrix,
                (
                    const lduMatrix& matrix,
                    const FieldField<Field, scalar>& coupleBouCoeffs,
                    const FieldField<Field, scalar>& coupleIntCoeffs,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& dict
                ),
                (
                    matrix,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaces,
                    dict
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                asymMatrix,
                (
                    const lduMatrix& matrix,
                    const FieldField<Field, scalar>& coupleBouCoeffs,
                    const FieldField<Field, scalar>& coupleIntCoeffs,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& dict
                ),
                (
                    matrix,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaces,
                    dict
                )
            );


        // Constructors

            preconditioner
            (
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces
            );


        // Selectors

            //- Return a new preconditioner
            static autoPtr<preconditioner> New
            (
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& dict
            );


        // Destructor

            virtual ~preconditioner()
            {}


        // Member functions

            //- Read and reset the preconditioner parameters
            //  from the given stream
            virtual void read(const dictionary&)
            {}

            //- Return wA the preconditioned form of residual rA
            virtual void precondition
            (
                scalarField& wA,
                const scalarField& rA,
                const direction cmpt = 0
            ) const = 0;

            //- Return wT the transpose-matrix preconditioned form of
            //  residual rT.  Used for asymmetric preconditioners
            virtual void preconditionT
            (
                scalarField& wT,
                const scalarField& rT,
                const direction cmpt = 0
            ) const
            {
                notImplemented
                (
                    type() +"::preconditionT"
                    "(scalarField& wT, const scalarField& rT, "
                    "const direction cmpt)"
                );
            }
    };


    // Static data

        // Declare name of the class and its debug switch
        ClassName("lduMatrix");

        //- Large scalar for the use in solvers
        static const scalar great_;

        //- Small scalar for the use in solvers
        static const scalar small_;


    // Constructors

        //- Construct given an LDU addressed mesh.
        //  The coefficients are initially empty for subsequent setting.
        explicit lduMatrix(const lduMesh&);

        //- Construct as copy
        lduMatrix(const lduMatrix&);

        //- Construct as copy or re-use as specified.
        lduMatrix(lduMatrix&, bool reUse);

        //- Construct given an LDU addressed mesh and an Istream
        //  from which the coefficients are read
        lduMatrix(const lduMesh&, Istream&);


        //- Clone
        autoPtr<lduMatrix> clone() const
        {
            return autoPtr<lduMatrix>
            (
                new lduMatrix(*this)
            );
        }


    // Destructor

        virtual ~lduMatrix();


    // Member functions

        // Access to addressing

            //- Return the LDU mesh from which the addressing is obtained
            const lduMesh& mesh() const
            {
                return lduMesh_;
            }

            //- Return the LDU addressing
            const lduAddressing& lduAddr() const
            {
                return lduMesh_.lduAddr();
            }

            //- Return the patch evaluation schedule
            const lduSchedule& patchSchedule() const
            {
                return lduAddr().patchSchedule();
            }


        // Access to coefficients

            scalarField& lower();
            scalarField& diag();
            scalarField& upper();

            const scalarField& lower() const;
            const scalarField& diag() const;
            const scalarField& upper() const;

            bool hasDiag() const
            {
                return (diagPtr_);
            }

            bool hasUpper() const
            {
                return (upperPtr_);
            }

            bool hasLower() const
            {
                return (lowerPtr_);
            }

            bool diagonal() const
            {
                return (diagPtr_ && !lowerPtr_ && !upperPtr_);
            }

            bool symmetric() const
            {
                return (diagPtr_ && (!lowerPtr_ && upperPtr_));
            }

            bool asymmetric() const
            {
                return (diagPtr_ && lowerPtr_ && upperPtr_);
            }


        // operations

            void sumDiag();
            void negSumDiag();

            void sumMagOffDiag(scalarField& sumOff) const;

            //- Matrix multiplication with updated interfaces.
            void Amul
            (
                scalarField&,
                const scalarField&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&,
                const direction cmpt
            ) const;

            //- Matrix multiplication without interfaces
            //  Result will be added to Ax
            void AmulCore
            (
                scalarField& Ax,
                const scalarField& x
            ) const;


            //- Matrix transpose multiplication with updated interfaces.
            void Tmul
            (
                scalarField&,
                const scalarField&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&,
                const direction cmpt
            ) const;

            //- Matrix transpose multiplication with updated coupled interfaces
            //  Result will be added to Tx
            void TmulCore
            (
                scalarField& Tx,
                const scalarField& x
            ) const;


            //- Sum the coefficients on each row of the matrix
            void sumA
            (
                scalarField&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&
            ) const;


            void residual
            (
                scalarField& rA,
                const scalarField& x,
                const scalarField& b,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const direction cmpt
            ) const;

            tmp<scalarField> residual
            (
                const scalarField& x,
                const scalarField& b,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const direction cmpt
            ) const;


            //- Initialise the update of coupled interfaces
            //  for matrix operations
            void initMatrixInterfaces
            (
                const FieldField<Field, scalar>& interfaceCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const scalarField& xif,
                scalarField& result,
                const direction cmpt
            ) const;

            //- Update coupled interfaces for matrix operations
            void updateMatrixInterfaces
            (
                const FieldField<Field, scalar>& interfaceCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const scalarField& xif,
                scalarField& result,
                const direction cmpt
            ) const;

            //- Initialise the update of coupled interfaces
            //  for matrix operations using buffered unscheduled transfer
            void bufferedInitMatrixInterfaces
            (
                const FieldField<Field, scalar>& interfaceCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const scalarField& xif,
                scalarField& result,
                const direction cmpt
            ) const;

            //- Update coupled interfaces for matrix operations
            //  using buffered unscheduled transfer
            void bufferedUpdateMatrixInterfaces
            (
                const FieldField<Field, scalar>& interfaceCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const scalarField& xif,
                scalarField& result,
                const direction cmpt
            ) const;


            template<class Type>
            tmp<Field<Type> > H(const Field<Type>&) const;

            template<class Type>
            tmp<Field<Type> > H(const tmp<Field<Type> >&) const;

            tmp<scalarField> H1() const;

            template<class Type>
            tmp<Field<Type> > faceH(const Field<Type>&) const;

            template<class Type>
            tmp<Field<Type> > faceH(const tmp<Field<Type> >&) const;


    // Member operators

        void operator=(const lduMatrix&);

        void negate();

        void operator+=(const lduMatrix&);
        void operator-=(const lduMatrix&);

        void operator*=(const scalarField&);
        void operator*=(scalar);


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const lduMatrix&);
};


// Helper typedefs
typedef lduMatrix::solver lduSolver;
typedef lduMatrix::solverPerformance lduSolverPerformance;
typedef lduMatrix::preconditioner lduPreconditioner;
typedef lduMatrix::smoother lduSmoother;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "lduMatrixTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
