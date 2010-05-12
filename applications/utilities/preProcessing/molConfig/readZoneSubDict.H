
// Info << "Zone description subDict " << cZ <<": "  << cellZoneI[cZ].name() << endl;

const dictionary& subDictI =
    molConfigDescription_.subDict(cellZoneI[cZ].name());

const scalar temperature(readScalar(subDictI.lookup("temperature")));

const word velocityDistribution(subDictI.lookup("velocityDistribution"));

const vector bulkVelocity(subDictI.lookup("bulkVelocity"));

const word id(subDictI.lookup("id"));

const scalar mass(readScalar(subDictI.lookup("mass")));

scalar numberDensity_read(0.0);

if (subDictI.found("numberDensity"))
{
    numberDensity_read = readScalar(subDictI.lookup("numberDensity"));
}
else if (subDictI.found("massDensity"))
{
    numberDensity_read = readScalar(subDictI.lookup("massDensity"))/mass;
}
else
{
    FatalErrorIn("readZoneSubDict.H\n")
        << "massDensity or numberDensity not specified " << nl
        << abort(FatalError);
}

const scalar numberDensity(numberDensity_read);

const word latticeStructure(subDictI.lookup("latticeStructure"));

const vector anchorPoint(subDictI.lookup("anchor"));

const word originSpecifies(subDictI.lookup("anchorSpecifies"));

if
(
    originSpecifies != "corner"
 && originSpecifies != "molecule"
)
{
    FatalErrorIn("readZoneSubDict.H\n")
        << "anchorSpecifies must be either 'corner' or 'molecule', found "
        << originSpecifies << nl
        << abort(FatalError);
}

bool tethered = false;

if (subDictI.found("tethered"))
{
    tethered = Switch(subDictI.lookup("tethered"));
}

const vector orientationAngles(subDictI.lookup("orientationAngles"));

scalar phi(orientationAngles.x()*mathematicalConstant::pi/180.0);
scalar theta(orientationAngles.y()*mathematicalConstant::pi/180.0);
scalar psi(orientationAngles.z()*mathematicalConstant::pi/180.0);

const tensor latticeToGlobal
(
    cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
    cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
    sin(psi)*sin(theta),
  - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
  - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
    cos(psi)*sin(theta),
    sin(theta)*sin(phi),
  - sin(theta)*cos(phi),
    cos(theta)
);

// Info << "\tcells: " << cellZoneI[cZ].size() << endl;
// Info << "\tnumberDensity: " << numberDensity << endl;
// Info << "\ttemperature: " << temperature << endl;
// Info << "\tvelocityDistribution: " << velocityDistribution << endl;
// Info << "\tbulkVelocity: " << bulkVelocity << endl;
// Info << "\tid: " << id << endl;
// Info << "\tmass: " << mass << endl;
// Info << "\tlatticeStructure: " << latticeStructure << endl;
// Info << "\tanchor: " << anchorPoint << endl;
// Info << "\toriginSpecifies: " << originSpecifies << endl;
// Info << "\ttethered: " << tethered << endl;
// Info << "\torientationAngles: " << orientationAngles << endl;
// Info << "\tlatticeToGlobal: " << latticeToGlobal << endl;

