    // Private data related to scalar functions

        //- Function pointer (used to avoid a conditional at every evaluation
        //  for debug levels > 0)
        void (Foam::equationReader::*reportEmbeddedDispatchFunction_)() const;

        //- Function pointer (used to avoid a conditional at every evaluation
        //  for debug levels > 0)
        void (Foam::equationReader::*reportEmbeddedReturnFunction_)() const;

        //- Function pointer (used to avoid a conditional at every evaluation
        //  for debug levels 1,2,5 or 6)
        void (Foam::equationReader::*reportScalarEvalStartFunction_)
        (
            const label&
        ) const;

        //- Function pointer (used to avoid a conditional at every operation
        //  for debug level 2 or 6)
        void (Foam::equationReader::*reportScalarOperationFunction_)
        (
            const label&,
            const label&
        ) const;

        //- Function pointer (used to avoid a conditional at every operation
        //  for debug level 2 or 6)
        void (Foam::equationReader::*reportScalarResultFunction_)
        (
            const scalar&
        ) const;

        //- Function pointer (used to avoid a conditional at every evaluation
        //  for debug levels 1,2,5 or 6)
        void (Foam::equationReader::*reportScalarEvalEndFunction_)
        (
            const scalar&
        ) const;

    // Private data related to dimensionSet functions

        //- Function pointer (used to avoid a conditional at every evaluation
        //  for debug levels 3,4,5 or 6)
        void (Foam::equationReader::*reportDimsEvalStartFunction_)
        (
            const label&
        ) const;

        //- Function pointer (used to avoid a conditional at every operation
        //  for debug level 4 or 6)
        void (Foam::equationReader::*reportDimsOperationFunction_)
        (
            const label&,
            const label&
        ) const;

        //- Function pointer (used to avoid a conditional at every operation
        //  for debug level 4 or 6)
        void (Foam::equationReader::*reportDimsResultFunction_)
        (
            const dimensionSet&
        ) const;

        //- Function pointer (used to avoid a conditional at every evaluation
        //  for debug levels 3,4,5 or 6)
        void (Foam::equationReader::*reportDimsEvalEndFunction_)
        (
            const dimensionSet&
        ) const;

    // Private member functions related to scalar functions

        //- Pointed to by reportEmbeddedDispatchFunction_ when debug level
        //  is 0 - this function does nothing
        void reportEmbeddedDispatchDisabled () const;

        //- Pointed to by reportEmbeddedDispatchFunction_ when debug level
        //  is not 0 - this function prints a line to the console
        void reportEmbeddedDispatchEnabled () const;

        //- Pointed to by reportEmbeddedReturnFunction_ when debug level
        //  is 0 - this function does nothing
        void reportEmbeddedReturnDisabled () const;

        //- Pointed to by reportEmbeddedReturnFunction_ when debug level
        //  is not 0 - this function prints a line to the console
        void reportEmbeddedReturnEnabled () const;

        //- Pointed to by reportScalarEvalStartFunction_ when debug level
        //  is none of: 1,2,5,6 - this function does nothing
        void reportScalarEvalStartDisabled
        (
            const label& index
        ) const;

        //- Pointed to by reportScalarEvalStartFunction_ when debug level
        //  is any of: 1,2,5,6 - this function reports evaluation has
        //  commenced (with details) to the console
        void reportScalarEvalStartEnabled
        (
            const label& index
        ) const;

        //- Pointed to by reportScalarOperationFunction_ when debug level
        //  is neither 2 nor 6 - this function does nothing
        void reportScalarOperationDisabled
        (
            const label& index,
            const label& i
        ) const;

        //- Pointed to by reportScalarOperationFunction_ when debug level
        //  is 2 or 6.  This function sends operation-by-operation
        //  information to the console
        void reportScalarOperationEnabled
        (
            const label& index,
            const label& i
        ) const;

        //- Pointed to by reportScalarResultFunction_ when debug level is
        //  neither 2 nor 6 - this function does nothing
        void reportScalarResultDisabled
        (
            const scalar& x
        ) const;

        //- Pointed to by reportScalarResultFunction_ when debug level is
        //  2 or 6.  This function reports the result operation-by-
        //  operation to the console
        void reportScalarResultEnabled
        (
            const scalar& x
        ) const;

        //- Pointed to by reportScalarEvalEndFunction_ when debug level
        //  is none of: 1,2,5,6 - this function does nothing
        void reportScalarEvalEndDisabled
        (
            const scalar& x
        ) const;

        //- Pointed to by reportScalarEvalEndFunction_ when debug level
        //  is any of: 1,2,5,6 - this function reports evaluation has
        //  completed (with details) to the console
        void reportScalarEvalEndEnabled
        (
            const scalar& x
        ) const;

    // Private member functions related to dimensionSet function

        //- Pointed to by reportDimsEvalStartFunction_ when debug level
        //  is none of: 3,4,5,6 - this function does nothing
        void reportDimsEvalStartDisabled
        (
            const label& index
        ) const;

        //- Pointed to by reportDimsEvalStartFunction_ when debug level
        //  is any of: 3,4,5,6 - this function reports evaluation has
        //  commenced (with details) to the console
        void reportDimsEvalStartEnabled
        (
            const label& index
        ) const;

        //- Pointed to by reportDimsOperationFunction_ when debug level
        //  is neither 4 nor 6 - this function does nothing
        void reportDimsOperationDisabled
        (
            const label& index,
            const label& i
        ) const;

        //- Pointed to by reportDimsOperationFunction_ when debug level
        //  is 4 or 6.  This function sends operation-by-operation
        //  information to the console
        void reportDimsOperationEnabled
        (
            const label& index,
            const label& i
        ) const;

        //- Pointed to by reportDimsResultFunction_ when debug level is
        //  neither 4 nor 6- this function does nothing
        void reportDimsResultDisabled
        (
            const dimensionSet& xDims
        ) const;

        //- Pointed to by reportDimsResultFunction_ when debug level is
        //  4 or 6.  This function reports the result operation-by-
        //  operation to the console
        void reportDimsResultEnabled
        (
            const dimensionSet& xDims
        ) const;

        //- Pointed to by reportDimsEvalEndFunction_ when debug level
        //  is none of: 3,4,5,6 - this function does nothing
        void reportDimsEvalEndDisabled
        (
            const dimensionSet& xDims
        ) const;

        //- Pointed to by reportDimsEvalEndFunction_ when debug level
        //  is any of: 3,4,5,6 - this function reports evaluation has
        //  completed (with details) to the console
        void reportDimsEvalEndEnabled
        (
            const dimensionSet& xDims
        ) const;

