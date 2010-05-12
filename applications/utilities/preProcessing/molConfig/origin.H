// Please refer to notes

// 1. Determine the unit cell dimensions: xU, yU and zU

const scalar xU = gap.x();
const scalar yU = gap.y();
const scalar zU = gap.z();

// 2. Determine the anchorPoint co-ordinates: xA, yA  and zA

const scalar xA = anchorPoint.x();
const scalar yA = anchorPoint.y();
const scalar zA = anchorPoint.z();

// 3. Determine the vector rAB from global co-ordinate system:

const vector rAB((xMid - xA), (yMid - yA), (zMid - zA));

// 4. Transform vector rAS into lattice co-ordinate system:

const vector rASTransf = transform(latticeToGlobal.T(), rAB);

// Info << "The vector rAS = " << rAS << endl;
// Info << "The vector rAStransf = " << rAStransf << endl;

// 5. Calculate the integer values: ni, nj and nk
scalar nIscalar = rASTransf.x()/xU;
scalar nJscalar = rASTransf.y()/yU;
scalar nKscalar = rASTransf.z()/zU;

// Info << "The nI, nJ, nK values before are: " << nIscalar <<" "<< nJscalar <<" "<< nKscalar << endl;

label nI = label(nIscalar + 0.5*sign(nIscalar));
label nJ = label(nJscalar + 0.5*sign(nJscalar));
label nK = label(nKscalar + 0.5*sign(nKscalar));

// Info << "The nI, nJ, nK values after are: " << nI <<" "<< nJ <<" "<< nK << endl;

// 6. Calculate the corrected starting point, rAC (in the lattice co-ordinate system):
const vector rAC((nI*xU), (nJ*yU), (nK*zU));

// 7. Transform the corrected starting point in the global co-ordinate system, rC:
const vector rC = anchorPoint + transform(latticeToGlobal, rAC);


const vector& origin = rC;

// Pout << "The Corrected Starting Point: " << origin << endl;

