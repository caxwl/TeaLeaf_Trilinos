#ifndef STEM_TRILINOS_H_
#define STEM_TRILINOS_H_

#include "Teuchos_RCPDecl.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_DefaultPlatform.hpp"

class TrilinosStem {
    typedef double Scalar;
    typedef int Ordinal;
    typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
    typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
    typedef Tpetra::Map<Ordinal,Ordinal,Node> Map;
    typedef Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal,Node> Matrix;
    typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node> Vector;
    typedef Tpetra::Operator<Scalar,Ordinal,Ordinal,Node> Operator;
    typedef Tpetra::MultiVector<Scalar,Ordinal,Ordinal,Node> MultiVector;

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
        static Teuchos::RCP<const Map> map;
        static Teuchos::RCP<Matrix> A;
        static Teuchos::RCP<Vector> b;
        static Teuchos::RCP<Vector> x;

        static Teuchos::RCP<Belos::LinearProblem<double, MultiVector, Operator> > problem;
        static Teuchos::RCP<Belos::SolverManager<double, MultiVector, Operator> > solver;

        static int* myGlobalIndices_;
        static int numLocalElements_;
};

#endif
