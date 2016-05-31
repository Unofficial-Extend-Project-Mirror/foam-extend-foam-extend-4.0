    sumLocalContErr =
        runTime.deltaT().value()*
        mag
        (
            fvc::ddt(rho)
          + fvc::div(phi)
        )().weightedAverage(rho*mesh.V()).value();

    globalContErr =
          runTime.deltaT().value()*
          (
              fvc::ddt(rho)
            + fvc::div(phi)
          )().weightedAverage(rho*mesh.V()).value();

    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;
