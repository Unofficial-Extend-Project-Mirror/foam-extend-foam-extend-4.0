#include "createTimeControls.H"

bool correctPhi
(
    piso.dict().lookupOrDefault("correctPhi", false)
);

bool checkMeshCourantNo
(
    piso.dict().lookupOrDefault("checkMeshCourantNo", false)
);
