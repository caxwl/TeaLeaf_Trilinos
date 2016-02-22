#ifndef TRILINOS_STEM_H_
#define TRILINOS_STEM_H_

#include "Teuchos_RCPDecl.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_DefaultPlatform.hpp"

#include "Ifpack2_Diagonal.hpp"

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
                int local_ymin, int local_ymax,
                int tl_max_iters, double tl_eps);
        static void solve(int nx, int ny,
                int local_xmin, int local_xmax,
                int local_ymin, int local_ymax,
                int global_xmin, int global_xmax,
                int global_ymin, int global_ymax,
                int* numIters,
                double rx,
                double ry,
                double* Kx,
                double* Ky,
                double* u0);
        static void finalise();
    private:
        static Teuchos::RCP<const Map> map;
        static Teuchos::RCP<Matrix> A;
        static Teuchos::RCP<Vector> b;
        static Teuchos::RCP<Vector> x;

        static Teuchos::RCP<Teuchos::ParameterList> solverParams;
        static Teuchos::RCP<Belos::LinearProblem<double, MultiVector, Operator> > problem;
        static Teuchos::RCP<Belos::SolverManager<double, MultiVector, Operator> > solver;
        static Teuchos::RCP<Ifpack2::Preconditioner<Scalar, Ordinal, Ordinal, Node> > preconditioner;

        static int* myGlobalIndices_;
        static int numLocalElements_;
        static int MyPID;
};

#endif
