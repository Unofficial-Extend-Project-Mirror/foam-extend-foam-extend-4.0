//======================================================================
//  variable 1 - var[0] , i.e variables are zero based
//======================================================================
int USERD_get_gold_variable_info
(
    char **var_description,
    char **var_filename,
    int *var_type,
    int *var_classify,
    int *var_complex,
    char **var_ifilename,
    float *var_freq,
    int *var_contran,
    int *var_timeset
)
{
#ifdef ENSIGHTDEBUG
    Info << "Entering: USERD_get_gold_variable_info" << endl 
        << flush;
#endif

    label offset = Num_variables - nSprayVariables 
        - nSecondMeshVariables - nFaMeshVariables;

    // scalars first ...
    for (label n=0; n<offset; n++)
    {
        if (isScalar[var2field[n]])
        {
            var_type[n]     = Z_SCALAR;
            var_classify[n] = Z_PER_ELEM;
            var_complex[n]  = FALSE;
            var_timeset[n]  = 1;
            strncpy
            (
                var_description[n],
                fieldNames[var2field[n]].c_str(),
                Z_BUFL
            );
        }
    }

    // ... and then vectors
    for (label n=0; n<offset; n++)
    {
        if (isVector[var2field[n]])
        {
            var_type[n]     = Z_VECTOR;
            var_classify[n] = Z_PER_ELEM;
            var_complex[n]  = FALSE;
            var_timeset[n]  = 1;
            strncpy
            (
                var_description[n],
                fieldNames[var2field[n]].c_str(),
                Z_BUFL
            );
        }
    }

    // ... and tensors (NB! all tensors are treated as asymmetric)
    for (label n=0; n<offset; n++)
    {
        if (isTensor[var2field[n]])
        {
            var_type[n]     = Z_TENSOR9;
            var_classify[n] = Z_PER_ELEM;
            var_complex[n]  = FALSE;
            var_timeset[n]  = 1;
            strncpy
            (
                var_description[n],
                fieldNames[var2field[n]].c_str(),
                Z_BUFL
            );
        }
    }

//     if (Numparts_available > nPatches+1)
    if (sprayPtr)
    {
        label Ns = lagrangianScalarNames.size();

        for (label n=0; n<Ns; n++)
        {
            var_type[offset + n] = Z_SCALAR;
            var_classify[offset + n] = Z_PER_ELEM;
            var_complex[offset + n] = FALSE;
            var_timeset[offset + n] = 1;
            word name = parcelPrepend + lagrangianScalarNames[n];
            strncpy
            (
                var_description[offset + n],
                name.c_str(),
                Z_BUFL
            );
        }

        for (label n=0; n<lagrangianVectorNames.size(); n++)
        {
            var_type[offset + Ns + n] = Z_VECTOR;
            var_classify[offset + Ns + n] = Z_PER_ELEM;
            var_complex[offset + Ns + n] = FALSE;
            var_timeset[offset + Ns  + n] = 1;
            word name = parcelPrepend + lagrangianVectorNames[n];
            strncpy
            (
                var_description[offset + Ns + n],
                name.c_str(),
                Z_BUFL
            );
        }
    }



    // Second mesh variables
    if(secondMeshPtr)
    {
    label start = Num_variables - nSecondMeshVariables - nFaMeshVariables;

    label offset = Num_variables - nFaMeshVariables;

    // scalars first ...
    for (label n=start; n<offset; n++)
    {
        if (secondMeshIsScalar[secondMeshVar2field[n-start]])
        {
            var_type[n]     = Z_SCALAR;
            var_classify[n] = Z_PER_ELEM;
            var_complex[n]  = FALSE;
            var_timeset[n]  = 1;
            strncpy
            (
                var_description[n],
                secondMeshFieldNames[secondMeshVar2field[n-start]].c_str(),
                Z_BUFL
            );
        }
    }

    // ... and then vectors
    for (label n=start; n<offset; n++)
    {
        if (secondMeshIsVector[secondMeshVar2field[n-start]])
        {
            var_type[n]     = Z_VECTOR;
            var_classify[n] = Z_PER_ELEM;
            var_complex[n]  = FALSE;
            var_timeset[n]  = 1;
            strncpy
            (
                var_description[n],
                secondMeshFieldNames[secondMeshVar2field[n-start]].c_str(),
                Z_BUFL
            );
        }
    }

    // ... and tensors (NB! all tensors are treated as asymmetric)
    for (label n=start; n<offset; n++)
    {
        if (secondMeshIsTensor[secondMeshVar2field[n-start]])
        {
            var_type[n]     = Z_TENSOR9;
            var_classify[n] = Z_PER_ELEM;
            var_complex[n]  = FALSE;
            var_timeset[n]  = 1;
            strncpy
            (
                var_description[n],
                secondMeshFieldNames[secondMeshVar2field[n-start]].c_str(),
                Z_BUFL
            );
        }
    }

    }


    // Finite area mesh variables
    if(faMeshPtr)
    {
        label start = Num_variables - nFaMeshVariables;

        label offset = Num_variables;

        // scalars first ...
        for (label n=start; n<offset; n++)
        {
            if (faMeshIsScalar[faMeshVar2field[n-start]])
            {
                var_type[n]     = Z_SCALAR;
                var_classify[n] = Z_PER_ELEM;
                var_complex[n]  = FALSE;
                var_timeset[n]  = 1;
                strncpy
                (
                    var_description[n],
                    faMeshFieldNames[faMeshVar2field[n-start]].c_str(),
                    Z_BUFL
                );

                Info << "Scalar: " << var_description[n] << endl;
            }
        }

        // ... and then vectors
        for (label n=start; n<offset; n++)
        {
            if (faMeshIsVector[faMeshVar2field[n-start]])
            {
                var_type[n]     = Z_VECTOR;
                var_classify[n] = Z_PER_ELEM;
                var_complex[n]  = FALSE;
                var_timeset[n]  = 1;
                strncpy
                (
                    var_description[n],
                    faMeshFieldNames[faMeshVar2field[n-start]].c_str(),
                    Z_BUFL
                );

                Info << "Vector: " << var_description[n] << endl;
            }
        }

        // ... and tensors (NB! all tensors are treated as asymmetric)
        for (label n=start; n<offset; n++)
        {
            if (faMeshIsTensor[faMeshVar2field[n-start]])
            {
                var_type[n]     = Z_TENSOR9;
                var_classify[n] = Z_PER_ELEM;
                var_complex[n]  = FALSE;
                var_timeset[n]  = 1;
                strncpy
                (
                    var_description[n],
                    faMeshFieldNames[faMeshVar2field[n-start]].c_str(),
                    Z_BUFL
                );

                Info << "Tensor: " << var_description[n] << endl;
            }
        }
    }

#ifdef ENSIGHTDEBUG
    Info << "Leaving: USERD_get_gold_variable_info" << endl 
        << flush;
#endif

    return Z_OK;
}



