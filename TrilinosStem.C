#include "TrilinosStem.h"

#include "BelosSolverFactory.hpp"
#include <BelosTpetraAdapter.hpp>

#include "Epetra_Map.h"

// Teuchos
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"

#define ARRAY2D(i,j,imin,jmin,ni) (i-(imin))+(((j)-(jmin))*(ni))

Teuchos::RCP<const TrilinosStem::Map> TrilinosStem::map;
Teuchos::RCP<TrilinosStem::Matrix> TrilinosStem::A;
Teuchos::RCP<TrilinosStem::Vector> TrilinosStem::b;
Teuchos::RCP<TrilinosStem::Vector> TrilinosStem::x;
int* TrilinosStem::myGlobalIndices_;
int TrilinosStem::numLocalElements_;


Teuchos::RCP<Belos::LinearProblem<double, TrilinosStem::MultiVector, TrilinosStem::Operator> > TrilinosStem::problem;
Teuchos::RCP<Belos::SolverManager<double, TrilinosStem::MultiVector, TrilinosStem::Operator> > TrilinosStem::solver;
Teuchos::RCP<Ifpack2::Preconditioner<TrilinosStem::Scalar, TrilinosStem::Ordinal, TrilinosStem::Ordinal, TrilinosStem::Node> > TrilinosStem::preconditioner;

extern "C" {
    void setup_trilinos_(
            int* nx,
            int* ny,
            int* localNx,
            int* localNy,
            int* local_xmin,
            int* local_xmax,
            int* local_ymin,
            int* local_ymax)
    {
        TrilinosStem::initialise(*nx, *ny,
                *localNx, *localNy,
                *local_xmin,
                *local_xmax,
                *local_ymin,
                *local_ymax);
    }

    void trilinos_solve_(
            int* nx,
            int* ny,
            int* local_xmin,
            int* local_xmax,
            int* local_ymin,
            int* local_ymax,
            int* global_xmin,
            int* global_xmax,
            int* global_ymin,
            int* global_ymax,
            double* rx,
            double* ry,
            double* Kx,
            double* Ky,
            double* u0)
    {
        TrilinosStem::solve(
                *nx,
                *ny,
                *local_xmin,
                *local_xmax,
                *local_ymin,
                *local_ymax,
                *global_xmin,
                *global_xmax,
                *global_ymin,
                *global_ymax,
                *rx,
                *ry,
                Kx,
                Ky,
                u0);
    }
}

void TrilinosStem::initialise(
        int nx,
        int ny,
        int localNx,
        int localNy,
        int local_xmin,
        int local_xmax,
        int local_ymin,
        int local_ymax)
{
    std::cout << "[STEM]: Setting up Trilinos/Tpetra...";

    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
    Teuchos::RCP<Node> node = platform.getNode();

    int numGlobalElements = nx * ny;
    numLocalElements_ = localNx * localNy;
    int rowSpacing = nx - localNx;

    myGlobalIndices_ = new int[numLocalElements_];

    int n = 0;

    int index = local_xmin;
    int i = 0;

    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            myGlobalIndices_[i] = index;
            i++;
            index++;
        }
        index += rowSpacing;
    }

    std::cerr << "\t[STEM]: creating map...";
    Teuchos::ArrayView<const Ordinal> localElementList = Teuchos::ArrayView<const Ordinal>(myGlobalIndices_, numLocalElements_);
    map = Teuchos::rcp(new Map(numGlobalElements, localElementList, 1, comm, node));
    std::cerr << " DONE. " << std::endl;


    size_t* numNonZero = new size_t[numLocalElements_];

    i = 0;

    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            size_t nnz = 1;

            if(1 != k)
                nnz++;
            if(ny != k)
                nnz++;
            if(1 != j)
                nnz++;
            if(nx != j)
                nnz++;

            numNonZero[i] = nnz;
            i++;
        }
    }

    std::cerr << "\t[STEM]: creating CrsMatrix...";
    Teuchos::ArrayRCP<const size_t> nnz = Teuchos::ArrayRCP<const size_t>(numNonZero, 0, numLocalElements_, false);

    A = Teuchos::rcp(new Matrix(map, nnz, Tpetra::StaticProfile));
    std::cerr << " DONE." << std::endl;

    b = Teuchos::rcp(new Vector(map));
    x = Teuchos::rcp(new Vector(map));

    problem = Teuchos::rcp(new Belos::LinearProblem<Scalar, MultiVector, Operator>(A, x, b));

    Belos::SolverFactory<Scalar, MultiVector, Operator> factory;

    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
    solverParams->set("Maximum Iterations", 1000);
     // TODO: Tolerance level causing issues at scale??
    solverParams->set("Convergence Tolerance", 1.0e-10); 

    solver = factory.create("RCG", solverParams);
    std::cout << "DONE." << std::endl;
}

void TrilinosStem::solve(
        int nx,
        int ny,
        int local_xmin,
        int local_xmax,
        int local_ymin,
        int local_ymax,
        int global_xmin,
        int global_xmax,
        int global_ymin,
        int global_ymax,
        double rx,
        double ry,
        double* Kx,
        double* Ky,
        double* u0)
{
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    bool insertValues = false;

    if(A->isFillComplete()) {
        A->resumeFill();
    } else {
        insertValues = true;
    }

    //A->setAllToScalar(0.0);

    std::vector<double> Values(4);
    std::vector<int> Indices(4);

    int numEntries = 0;
    int local_nx = local_xmax - local_xmin + 1 + 4;

    int i = 0;

    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            Values.clear();
            Indices.clear();

            double c2 = -rx*Kx[ARRAY2D(j,k,local_xmin-2,local_ymin-2,local_nx)];
            double c3 = -rx*Kx[ARRAY2D(j+1,k,local_xmin-2,local_ymin-2,local_nx)];
            double c4 = -ry*Ky[ARRAY2D(j,k,local_xmin-2,local_ymin-2,local_nx)];
            double c5 = -ry*Ky[ARRAY2D(j,k+1,local_xmin-2,local_ymin-2,local_nx)];

            numEntries = 0;

            if(j == global_xmin) {
                c2 = 0;
            }
            if(1 != j) {
                numEntries++;
                Values.push_back(c2);
                Indices.push_back(myGlobalIndices_[i]-1);
            }


            if(j == global_xmax) {
                c3 = 0;
            }
            if(nx != j) {
                numEntries++;
                Values.push_back(c3);
                Indices.push_back(myGlobalIndices_[i]+1);
            }

            if (k == global_ymin) {
                c4 = 0;
            }
            if(1 != k) {
                numEntries++;
                Values.push_back(c4);
                Indices.push_back(myGlobalIndices_[i]-nx);
            }

            if (k == global_ymax) {
                c5 = 0;
            }
            if(ny != k) {
                numEntries++;
                Values.push_back(c5);
                Indices.push_back(myGlobalIndices_[i]+nx);
            }


            double diagonal = (1.0 -c2 -c3 -c4 -c5);


            if (insertValues) {
                A->insertGlobalValues(myGlobalIndices_[i],
                        Teuchos::ArrayView<Ordinal>(&Indices[0], numEntries), 
                        Teuchos::ArrayView<Scalar>(&Values[0], numEntries));

                A->insertGlobalValues(myGlobalIndices_[i],
                        Teuchos::tuple<Ordinal>( myGlobalIndices_[i] ),
                        Teuchos::tuple<Scalar>(diagonal));
            } else {
                A->replaceGlobalValues(myGlobalIndices_[i],
                        Teuchos::ArrayView<Ordinal>(&Indices[0], numEntries),
                        Teuchos::ArrayView<Scalar>(&Values[0], numEntries));

                A->replaceGlobalValues(myGlobalIndices_[i],
                        Teuchos::tuple<Ordinal>( myGlobalIndices_[i] ),
                        Teuchos::tuple<Scalar>(diagonal));
            }

            i++;
        }
    }

    A->fillComplete();

    Indices.clear();
    Values.clear();

    i = 0;
    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            double c2 = Kx[ARRAY2D(j,k,local_xmin-2,local_ymin-2,local_nx)];
            double c3 = Kx[ARRAY2D(j+1,k,local_xmin-2,local_ymin-2,local_nx)];
            double c4 = Ky[ARRAY2D(j,k,local_xmin-2,local_ymin-2,local_nx)];
            double c5 = Ky[ARRAY2D(j,k+1,local_xmin-2,local_ymin-2,local_nx)];

            double value = u0[ARRAY2D(j,k,local_xmin-2, local_ymin-2, local_nx)];

            b->replaceGlobalValue(myGlobalIndices_[i], value);
            i++;
        }
    }

    i = 0;
    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            double value = u0[ARRAY2D(j,k, local_xmin-2, local_ymin-2, local_nx)];

            x->replaceGlobalValue(myGlobalIndices_[i], value);
            i++;
        }
    }

    problem->setOperator(A);
    problem->setLHS(x);
    problem->setRHS(b);

    problem->setProblem();

    const Teuchos::RCP<const TrilinosStem::Matrix> const_ptr_to_A = A;
    preconditioner = Teuchos::rcp(new Ifpack2::Diagonal<const TrilinosStem::Matrix>(const_ptr_to_A));
    preconditioner->compute();
    problem->setLeftPrec(preconditioner);

    solver->setProblem(problem);

    Belos::ReturnType result = solver->solve();

    const int numIters = solver->getNumIters();
    std::cout << "[STEM]: num_iters = " << numIters << std::endl;

    Teuchos::Array<Scalar> solution(numLocalElements_);
    x->get1dCopy(solution, numLocalElements_);

    i = 0;
    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            u0[ARRAY2D(j,k,local_xmin-2, local_ymin-2, local_nx)] = solution[i];
            i++;
        }
    }
}
