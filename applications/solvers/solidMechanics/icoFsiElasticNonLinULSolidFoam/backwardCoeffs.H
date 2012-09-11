
dimensionedScalar Cn("Cn", dimless/dimTime, 0.0);
dimensionedScalar Co("Co", dimless/dimTime, 0.0);
dimensionedScalar Coo("Coo", dimless/dimTime, 0.0);

scalar deltaT = runTime.deltaT().value();
scalar deltaT0 = runTime.deltaT0().value();

Cn.value() = 1 + deltaT/(deltaT + deltaT0);
Coo.value() = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
Co.value() = Cn.value() + Coo.value();

if(runTime.timeIndex() == 1)
{
    Cn.value() = 1.0;
    Co.value() = 1.0;
    Coo.value() = 0.0;
}

Cn.value() /= deltaT;
Co.value() /= deltaT;
Coo.value() /= deltaT;
