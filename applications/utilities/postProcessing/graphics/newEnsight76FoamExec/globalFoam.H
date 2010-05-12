//======================================================================
// Global variables
const int maxNames = 1000;

// define the name to be displayed in the window.
static char readerName[]           = "newOpenFOAM";
static char meshName[]             = "cells";
static char readerVersion[]        = "2.03";

// everything is one part in foam, except the spray
static int Num_unstructured_parts  = 1;
static int Num_structured_parts    = 0;
static int Numparts_available      = 1;
static int nPatches                = 0;

static int Num_timesets            = 1;
static int Geom_timeset_number     = 1;
static int Num_time_steps          = 1;
// static int Num_global_nodes        = 0;
static int Num_variables           = 0;
static int Num_dataset_files       = 0;
static int Current_time_step       = 0;

static label nSprayVariables       = 0;
static label nMaxParcels           = 0;

static bool isScalar[maxNames];
static bool isVector[maxNames];
static bool isTensor[maxNames];
static bool isSpray[maxNames];

static word scalarName             = "volScalarField";
static word vectorName             = "volVectorField";
static word tensorName             = "volTensorField";
static word sprayName              = "lagrangian";
static word sprayScalarFieldName   = "scalarField";
static word sprayVectorFieldName   = "vectorField";
static word sprayTensorFieldName   = "tensorField";
static word parcelPrepend          = "parcel_";
static word pointPrepend           = "point_";

static fileName rootDir;
static fileName caseDir;

static instantList TimeList;

static List<word> fieldNames;
static List<word> lagrangianScalarNames;
static List<word> lagrangianVectorNames;
static label var2field[maxNames];

static Time *runTimePtr = 0;
static fvMesh *meshPtr = 0;
static Cloud<passiveParticle> *sprayPtr = 0;

// Second mesh global data
static fvMesh* secondMeshPtr = 0;
static int nSecondMeshPatches = 0;
static char secondMeshName[] = "cells_2";
static List<word> secondMeshFieldNames;
static label secondMeshVar2field[maxNames];
static bool secondMeshIsScalar[maxNames];
static bool secondMeshIsVector[maxNames];
static bool secondMeshIsTensor[maxNames];
static label nSecondMeshVariables = 0;


// Finite area mesh global data
static faMesh* faMeshPtr = 0;
static word faScalarName             = "areaScalarField";
static word faVectorName             = "areaVectorField";
static word faTensorName             = "areaTensorField";
static char faMeshName[] = "faCells";
static List<word> faMeshFieldNames;
static label faMeshVar2field[maxNames];
static bool faMeshIsScalar[maxNames];
static bool faMeshIsVector[maxNames];
static bool faMeshIsTensor[maxNames];
static label nFaMeshVariables = 0;
static label secondMeshPartNum = 0;


// Additional data
static float scalingFactor = 1.0;
