Pav << runTime.theta() << "   " << p.weightedAverage(mesh.V()).value() << endl;
Tav << runTime.theta() << "   " << T.weightedAverage(mesh.V()).value() << endl;

Info << "Max T = " << max(T) << " K" << endl;
Info << "Min T = " << min(T) << " K" << endl;
Info << "Max p = " << max(p)/1.0e5 << " bar" << endl;

massBal << runTime.theta() << "   "
        << fvc::domainIntegrate(rho).value() << endl;
volume << runTime.theta() << "   " << sum(mesh.V()) << endl;

debugT << nl << "Crank angle: " << runTime.theta() << endl;
debugT << "Max T = " << max(T) << " K" << endl;
debugT << "Min T = " << min(T) << " K" << endl;

kav << runTime.theta() << "   "
    << (turbulence->k())().weightedAverage(mesh.V()).value() << endl;

epsav << runTime.theta() << "   "
      << (turbulence->epsilon())().weightedAverage(mesh.V()).value() << endl;

debugMinT << runTime.theta() << "   " << min(T).value() << endl;
debugMaxT << runTime.theta() << "   " << max(T).value() << endl;
