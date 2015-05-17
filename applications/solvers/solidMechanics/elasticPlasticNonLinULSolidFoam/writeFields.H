if (runTime.outputTime())
{
    volScalarField epsilonEq
    (
        IOobject
        (
            "epsilonEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((2.0/3.0)*magSqr(dev(epsilon)))
    );

    Info<< "Max epsilonEq = " << max(epsilonEq).value() << endl;

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );

    Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;

    pointMesh pMesh(mesh);
    pointScalarField contactPointGap
    (
        IOobject
        (
            "contactPointGap",
            runTime.timeName(),
            mesh,
           IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(mesh.boundary(), patchi)
    {
        if(DU.boundaryField()[patchi].type() == solidContactFvPatchVectorField::typeName)
        {
            const solidContactFvPatchVectorField& DUpatch =
                refCast<const solidContactFvPatchVectorField>
                    (DU.boundaryField()[patchi]);

            if(!DUpatch.master())
            {
                const labelList& meshPoints = mesh.boundaryMesh()[patchi].meshPoints();
                const scalarField gap = DUpatch.normalContactModelPtr()->slaveContactPointGap();
                forAll(meshPoints, pointi)
                {
                    contactPointGap[meshPoints[pointi]] = gap[pointi];
                }
            }
        }
    }


    //- total force
/*
    forAll(mesh.boundary(), patchi)
    {
        Info << "Patch " << mesh.boundary()[patchi].name() << endl;
        vectorField totalForce = mesh.Sf().boundaryField()[patchi] & sigma.boundaryField()[patchi];

        vector force = sum( totalForce );
        Info << "\ttotal force is " << force << " N" << endl;

        tensorField F = I + gradDU.boundaryField()[patchi];
        tensorField Finv = inv(F);
        //vectorField nCurrent = Finv & n.boundaryField()[patchi];
        //nCurrent /= mag(nCurrent);
        //scalar normalForce = sum( nCurrent & totalForce );
        scalar normalForce = sum( n.boundaryField()[patchi] & totalForce );
        Info << "\tnormal force is " << normalForce << " N" << endl;
        //scalar shearForce = mag(sum( (I - sqr(nCurrent)) & totalForce ));
        scalar shearForce = mag(sum( (I - sqr(n.boundaryField()[patchi])) & totalForce ));
        Info << "\tshear force is " << shearForce << " N" << endl;

        // if(mesh.boundary()[patchi].type() != "empty")
        //   {
    //     vector Sf0 = Sf.boundaryField()[patchi][0];
    //     symmTensor sigma0 = sigma.boundaryField()[patchi][0];
    //     Info << "sigmab[0] is " << sigma0 << nl
    //      << "Sfb is " << Sf0 << nl
    //      << "force is " << (Sf.boundaryField()[patchi][0]&sigma.boundaryField()[patchi][0]) << nl
    //      << "Sfx*sigmaxx " << (Sf0[vector::X]*sigma0[symmTensor::XX]) << nl
    //      << "Sfy*sigmaxy " << (Sf0[vector::Y]*sigma0[symmTensor::XY]) << nl
    //      << "Sfx*sigmayx " << (Sf0[vector::X]*sigma0[symmTensor::XY]) << nl
    //      << "Sfy*sigmayy " << (Sf0[vector::Y]*sigma0[symmTensor::YY]) << nl
    //      << endl;
        //vector SfTL(-0.000137451, 0.00383599, -4.76878e-20);
        // vector SfTL = Finv[0] & vector(0,0.004,0);
        // Info << "SfTLx*sigmaxx " << (SfTL[vector::X]*sigma0[symmTensor::XX]) << nl
        //      << "SfTLy*sigmaxy " << (SfTL[vector::Y]*sigma0[symmTensor::XY]) << nl
        //      << "SfTLx*sigmayx " << (SfTL[vector::X]*sigma0[symmTensor::XY]) << nl
        //      << "SfTLy*sigmayy " << (SfTL[vector::Y]*sigma0[symmTensor::YY]) << nl
        //      << endl;
    //      }
    }
*/

    runTime.write();
  }
