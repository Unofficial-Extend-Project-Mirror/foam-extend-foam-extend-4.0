/*
freesteam - IAPWS-IF97 steam tables library
Copyright (C) 2004-2009  John Pye

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#define FREESTEAM_BUILDING_LIB
#include "zeroin.h"

#include <math.h>
#include <stdio.h>

#ifndef DBL_EPSILON
	#define DBL_EPSILON 2e-16
#endif

char zeroin_solve(ZeroInSubjectFunction *func, void *user_data, double lowerbound, double upperbound, double tol, double *solution, double *error){

    double a, b, c;	///<  Abscissae, descr. see above.
    double fa;      ///<  f(a)
    double fb;      ///<  f(b)
    double fc;      ///<  f(c)

    a = lowerbound;
    b = upperbound;
    fa = (*func)(a,user_data);
	    fb = (*func)(b,user_data);
    c = a;
    fc = fa;

	if(fa == 0.){
		*error = 0.; // used by getError
		*solution = a;
		//fprintf(stderr,"perfect solution\n");
		return 0;
	}

    //  Main iteration loop

    for (;;) {
	    double prev_step = b - a;    ///<  Distance from the last but one to the last approximation
	    double tol_act;              ///<  Actual tolerance
	    double p;                    ///<  Interpolation step is calculated in the form p/q; division
	    double q;                    ///<  operations is delayed until the last moment
	    double new_step;             ///<  Step at this iteration

	    if (fabs(fc) < fabs(fb)) {
		    a = b;
		    b = c;
		    c = a;  //  Swap data for b to be the best approximation
		    fa = fb;
		    fb = fc;
		    fc = fa;
	    }

	    // DBL_EPSILON is defined in math.h
	    tol_act = 2.0* DBL_EPSILON * fabs(b) + tol / 2.0;

	    new_step = (c - b) / 2.0;
		//fprintf(stderr,"step = %g\n",new_step);

	    if (fabs(new_step) <= tol_act || fb == 0.) {
			*error = fb;
			*solution = b;
		    //fprintf(stderr,"best solution is b: f(b=%g) = %g, f(a=%g) = %g\n",b,fb,a,fb);
		    return 0;
	    }
	    //  Decide if the interpolation can be tried

	    if (fabs(prev_step) >= tol_act   // If prev_step was large enough and was in true direction,
	        && fabs(fa) > fabs(fb))      // Interpolatiom may be tried
	    {
		    register double t1, t2;
		    double cb;

		    cb = c - b;
		    if (a == c) {
			    // If we have only two distinct points
			    // then only linear interpolation can be applied
			    t1 = fb / fa;
			    p = cb * t1;
			    q = 1.0 - t1;
		    } else {
			    // Quadric inverse interpolation

			    q = fa / fc;
			    t1 = fb / fc;
			    t2 = fb / fa;
			    p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
			    q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
		    }
		    if (p > 0.) {
			    // p was calculated with the opposite sign; make p positive-
			    q = -q;	// and assign possible minus to q
		    } else {
			    p = -p;
		    }

		    if (p < (0.75 * cb * q - fabs(tol_act * q) / 2.0)
		        && p < fabs(prev_step * q / 2.0)
		       ) {
			    // If b+p/q falls in [b,c] and
			    // isn't too large it is accepted
			    new_step = p / q;
		    }
		    // If p/q is too large then the bissection procedure can
		    // reduce [b,c] range to more extent
	    }

	    if (fabs(new_step) < tol_act) {	// Adjust the step to be not less
		    if (new_step > 0.)	// than tolerance
			    new_step = tol_act;
		    else
			    new_step = -tol_act;
	    }

	    a = b;
	    fa = fb;		// Save the previous approx.
	    b += new_step;
	    fb = (*func)(b,user_data);	//  Do step to a new approxim.

	    if ((fb > 0. && fc > 0.)
	        || (fb < 0. && fc < 0.)) {
		    c = a;
		    fc = fa;	//  Adjust c for it to have a sign opposite to that of b
	    }
    }
    // (((we never arrive here)))
}
