#ifndef STEM_TRILINOS_H_
#define STEM_TRILINOS_H_

#include "Teuchos_RCPDecl.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosSolverManager.hpp"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

class TrilinosStem {
    public:
        TrilinosStem();
        static void initialise( int nx, int ny,
                int localNx, int localNy,
                int local_xmin, int local_xmax,
                int local_ymin, int local_ymax);
        static void solve(int nx, int ny,
                int local_xmin, int local_xmax,
                int local_ymin, int local_ymax,
                int global_xmin, int global_xmax,
                int global_ymin, int global_ymax,
                double rx,
                double ry,
                double* Kx,
                double* Ky,
                double* u0);
    private:
        static Teuchos::RCP<Epetra_Map> map;
        static Teuchos::RCP<Epetra_CrsMatrix> A;
        static Teuchos::RCP<Epetra_Vector> b;
        static Teuchos::RCP<Epetra_Vector> x;

        static Teuchos::RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> > problem;
        static Teuchos::RCP<Belos::SolverManager<double, Epetra_MultiVector, Epetra_Operator> > solver;

        static int* myGlobalIndices_;
        static int numLocalElements_;
};

#endif
