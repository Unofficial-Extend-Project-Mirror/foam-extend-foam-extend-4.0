#include "point2D.H"
#include "HormannAgathos.H"

using namespace Foam;

// Self check ... for debugging
//
//
//                          w = 0
//
//
//   (-10, 10 X--------------X--------------X (10, 10)
//            |              .              |
//            |              .              |
//            |     w = 1    .  w = 1       |
//            |              .              |
//            |              .              |
//    w = 0   X..............X.(0,0)........X     w = 0
//            |              .              |
//            |              .              |
//            |     w = 1    .  w = 1       |
//            |              .              |
//            |              .              |
// (-10, -10) X--------------X--------------X (10, -10)
//
//
//                        w = 0
//
//
//


int main()
{
    Info << "Starting HormannAgathos::selfCheck" << endl;

    List<point2D> P(4);
    point2D p0(-10.0, -10);
    point2D p1( 10.0, -10);
    point2D p2( 10.0,  10);
    point2D p3(-10.0,  10);

    point2D deltaX = (p1 - p0);
    point2D deltaY = (p2 - p1);

    P[0] = p0;
    P[1] = p1;
    P[2] = p2;
    P[3] = p3;

    scalar distTol = 0.0000001;
    HormannAgathos ha(P, distTol);

    // Check every vertex
    ha.evaluate(p0);
    ha.evaluate(p1);
    ha.evaluate(p2);
    ha.evaluate(p3);

    // Ha.Evaluate on edges
    ha.evaluate(p0 + 0.5*deltaX);
    ha.evaluate(p1 + 0.5*deltaY);
    // In this case, we exit with ON_EDGE while w = 1; this is ok, we
    // chose to exit sooner
    ha.evaluate(p2 - 0.5*deltaX);
    ha.evaluate(p3 - 0.5*deltaY);
     
    // Ha.Evaluate right in the middle of each quadrant
    ha.evaluate(-0.25*deltaX - 0.25*deltaY);
    ha.evaluate( 0.25*deltaX - 0.25*deltaY);
    ha.evaluate( 0.25*deltaX + 0.25*deltaY);
    ha.evaluate(-0.25*deltaX + 0.25*deltaY);

    // Ha.Evaluate outside of polygon
    point2D hugeX(10000.0, 0.0);
    point2D hugeY( 0.0, 10000.0);
    ha.evaluate(p0 - hugeX - hugeY);
    ha.evaluate(p1 + hugeX - hugeY);
    ha.evaluate(p2 + hugeX + hugeY);
    ha.evaluate(p3 - hugeX + hugeY);
   
    ha.evaluate(p0 + 0.5*deltaX - hugeY);
    ha.evaluate(p1 + 0.5*deltaY + hugeX);
    ha.evaluate(p2 - 0.5*deltaX + hugeY);
    ha.evaluate(p3 - 0.5*deltaY - hugeX);

    // In problem:
    ha.evaluate(p3 - 0.5*deltaY);
    ha.evaluate(p3 - 0.5*deltaY - hugeX);
    ha.evaluate(p3 - 0.5*deltaY);
    ha.evaluate(p3 - 0.5*deltaY - hugeX);

    Info << "Exiting HormannAgathos selfCheck" << endl;
    return 0;
}
