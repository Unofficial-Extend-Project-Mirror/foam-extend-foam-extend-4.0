/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
   domain Java Matrix library JAMA. */

#include <math.h>
#include "fvCFD.H"

#include "eig3.H"

#include <Eigen/Dense>
using namespace Eigen;

#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))

#define n 3

static double hypot2(double x, double y) {
    return ::sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(tensor& V, vector& d, vector& e)
{
//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

    for (int j = 0; j < n; j++)
    {
        d[j] = V[6 + j];
    }

    // Householder reduction to tridiagonal form.

    for (int i = n-1; i > 0; i--)
    {
        // Scale to avoid under/overflow.

        double scale = 0.0;
        double h = 0.0;
        for (int k = 0; k < i; k++)
        {
            scale = scale + fabs(d[k]);
        }

        if (scale < SMALL)
//         if (scale == 0.0)
        {
            e[i] = d[i-1];
            for (int j = 0; j < i; j++)
            {
                d[j] = V[(i-1)*3 + j];
                V[i*3 + j] = 0.0;
                V[j*3 + i] = 0.0;
            }
        }
        else
        {
            // Generate Householder vector.

            for (int k = 0; k < i; k++)
            {
                d[k] /= scale;
                h += d[k] * d[k];
            }

            double f = d[i-1];
            double g = ::sqrt(h);
            if (f > 0)
            {
                g = -g;
            }

            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (int j = 0; j < i; j++)
            {
                e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.

            for (int j = 0; j < i; j++)
            {
                f = d[j];
                V[j*3 + i] = f;
                g = e[j] + V[j*3 + j] * f;
                for (int k = j+1; k <= i-1; k++)
                {
                    g += V[k*3 + j] * d[k];
                    e[k] += V[k*3 + j] * f;
                }
                e[j] = g;
            }

            f = 0.0;
            for (int j = 0; j < i; j++)
            {
                e[j] /= h;
                f += e[j] * d[j];
            }

            double hh = f / (h + h);
            for (int j = 0; j < i; j++)
            {
                e[j] -= hh * d[j];
            }

            for (int j = 0; j < i; j++)
            {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i-1; k++)
                {
                    V[k*3 + j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[(i-1)*3 + j];
                V[i*3 + j] = 0.0;
            }
        }

        d[i] = h;
    }

    // Accumulate transformations.

    for (int i = 0; i < n-1; i++)
    {
        V[(n-1)*3 + i] = V[i*3 + i];
        V[i*3 + i] = 1.0;
        double h = d[i+1];
        if (h != 0.0)
        {
            for (int k = 0; k <= i; k++)
            {
                d[k] = V[k*3 + i+1] / h;
            }
            for (int j = 0; j <= i; j++)
            {
                double g = 0.0;
                for (int k = 0; k <= i; k++)
                {
                    g += V[k*3 + i+1] * V[k*3 + j];
                }
                for (int k = 0; k <= i; k++)
                {
                    V[k*3 + j] -= g * d[k];
                }
            }
        }
        for (int k = 0; k <= i; k++)
        {
            V[k*3 + i+1] = 0.0;
        }
    }
    for (int j = 0; j < n; j++)
    {
        d[j] = V[(n-1)*3 + j];
        V[(n-1)*3 + j] = 0.0;
    }

    V[(n-1)*3 + n-1] = 1.0;
    e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.

static void tql2(tensor& V, vector& d, vector& e)
{
//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

    for (int i = 1; i < n; i++)
    {
        e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = ::pow(2.0,-52.0);
    for (int l = 0; l < n; l++)
    {
        // Find small subdiagonal element

        tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
        int m = l;
        while (m < n)
        {
            if (fabs(e[m]) <= eps*tst1)
            {
                break;
            }
            m++;
        }

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.

        if (m > l)
        {
            int iter = 0;

            do
            {
                iter = iter + 1;  // (Could check iteration count here.)

                // Compute implicit shift

                double g = d[l];
                double p = (d[l+1] - g) / (2.0 * e[l]);
                double r = hypot2(p,1.0);
                if (p < 0)
                {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                double dl1 = d[l+1];
                double h = g - d[l];
                for (int i = l+2; i < n; i++)
                {
                    d[i] -= h;
                }
                f = f + h;

                // Implicit QL transformation.

                p = d[m];
                double c = 1.0;
                double c2 = c;
                double c3 = c;
                double el1 = e[l+1];
                double s = 0.0;
                double s2 = 0.0;

                for (int i = m-1; i >= l; i--)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p,e[i]);
                    e[i+1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i+1] = h + s * (c * g + s * d[i]);

                    // Accumulate transformation.

                    for (int k = 0; k < n; k++)
                    {
                        h = V[k*3 + i+1];
                        V[k*3 + i+1] = s * V[k*3 + i] + c * h;
                        V[k*3 + i] = c * V[k*3 + i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                // Check for convergence.

            }
            while (fabs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.

    for (int i = 0; i < n-1; i++)
    {
        int k = i;
        double p = d[i];
        for (int j = i+1; j < n; j++)
        {
            if (d[j] < p)
            {
                k = j;
                p = d[j];
            }
        }
        if (k != i)
        {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++)
            {
                p = V[j*3 + i];
                V[j*3 + i] = V[j*3 + k];
                V[j*3 + k] = p;
            }
        }
    }
}


void eigen_decomposition(const tensor& A, tensor& V, vector& d)
{
    vector e = vector::zero;
    V = A;

    tred2(V, d, e);
    tql2(V, d, e);

    // transpose eigen vector
    V = V.T();
}


void eigen_decomposition(const symmTensor& sA, tensor& V, vector& d)
{
    vector e = vector::zero;
    V = tensor(sA);

    tred2(V, d, e);
    tql2(V, d, e);

    // transpose eigen vector
    V = V.T();



//     Info << "old " << d << endl;
//     Info << V.xx() << " " << V.xy() << " " << V.xz() << endl;
//     Info << V.yx() << " " << V.yy() << " " << V.yz() << endl;
//     Info << V.zx() << " " << V.zy() << " " << V.zz() << endl << endl;





//     Matrix3f A;

//     A << sA.xx(), sA.xy(), sA.xz(),
//         sA.xy(), sA.yy(), sA.yz(),
//         sA.xz(), sA.yz(), sA.zz();

//     SelfAdjointEigenSolver<Matrix3f> eigensolver(A);

//     if (eigensolver.info() != Success) abort();

//     V.xx() = eigensolver.eigenvectors()(0,0);
//     V.yx() = eigensolver.eigenvectors()(1,0);
//     V.zx() = eigensolver.eigenvectors()(2,0);

//     V.xy() = eigensolver.eigenvectors()(0,1);
//     V.yy() = eigensolver.eigenvectors()(1,1);
//     V.zy() = eigensolver.eigenvectors()(2,1);

//     V.xz() = eigensolver.eigenvectors()(0,2);
//     V.yz() = eigensolver.eigenvectors()(1,2);
//     V.zz() = eigensolver.eigenvectors()(2,2);

//     V = V.T();





//     std::cout << "\n" << A << "\n" << "\n";
//     std::cout << eigensolver.eigenvalues() << "\n" << "\n";
//     std::cout << eigensolver.eigenvectors() << "\n" << "\n";
}


void eigen_decomposition(const symmTensor& sA, tensor& V, diagTensor& d)
{
    vector e = vector::zero;
    V = tensor(sA);

    // We will use a temporary vector
    vector dVec = vector::zero;
    tred2(V, dVec, e);
    tql2(V, dVec, e);

    // And then transfer eigenvalues to diagTensor
    d.xx() = dVec.x();
    d.yy() = dVec.y();
    d.zz() = dVec.z();

    // transpose eigen vector
    V = V.T();
}
