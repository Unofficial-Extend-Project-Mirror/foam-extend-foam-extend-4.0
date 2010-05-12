        /*
        volTensorField gradU = fvc::grad(U);
        volSymmTensorField D = symm(fvc::grad(U));
        volTensorField Dprim = symm(fvc::grad(U - Umean));

        volScalarField prod = -((U - Umean)*(U - Umean)) && D;
        */

        /*
        volScalarField txx
        (
            IOobject
            (
                "txx",
                Times[i].name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 1, -1, 0, 0)
        );
        txx =sqrt(Txx - (Umeanx*Umeanx));
        txx.write();

        volScalarField tyy
        (
            IOobject
            (
                "tyy",
                Times[i].name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 1, -1, 0, 0)
        );
        tyy = sqrt(Tyy - (Umeany*Umeany));
        tyy.write();

        volScalarField tzz
        (
            IOobject
            (
                "tzz",
                Times[i].name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 1, -1, 0, 0)
        );
        tzz = sqrt(Tzz - (Umeanz*Umeanz));
        tzz.write();

        volScalarField txy
        (
            IOobject
            (
                "txy",
                Times[i].name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, 0, 0)
        );
        txy = Txy - (Umeanx*Umeany);
        txy.write();
*/
