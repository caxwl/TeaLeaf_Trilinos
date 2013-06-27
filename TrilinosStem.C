#include "TrilinosStem.h"

#include "BelosSolverFactory.hpp"
#include "Epetra_MpiComm.h"

#include "Teuchos_ParameterList.hpp"

#define ARRAY2D(i,j,imin,jmin,ni) (i-(imin))+(((j)-(jmin))*(ni))

Teuchos::RCP<Epetra_Map> TrilinosStem::map;
Teuchos::RCP<Epetra_CrsMatrix> TrilinosStem::A;
Teuchos::RCP<Epetra_Vector> TrilinosStem::b;
Teuchos::RCP<Epetra_Vector> TrilinosStem::x;
int* TrilinosStem::myGlobalIndices_;
int TrilinosStem::numLocalElements_;


Teuchos::RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> > TrilinosStem::problem;
Teuchos::RCP<Belos::SolverManager<double, Epetra_MultiVector, Epetra_Operator> > TrilinosStem::solver;

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
    std::cout << "[STEM]: Setting up Trilinos/Epetra...";

    Epetra_MpiComm comm(MPI_COMM_WORLD);

    long long numGlobalElements = nx * ny;
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
    map = Teuchos::rcp(new Epetra_Map(numGlobalElements, numLocalElements_, myGlobalIndices_, 1, comm));
    std::cerr << " DONE. " << std::endl;


    int* numNonZero = new int[numLocalElements_];

    i = 0;
    
    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            int nnz = 1;

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
    A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *map, numNonZero));
    std::cerr << " DONE." << std::endl;

    b = Teuchos::rcp(new Epetra_Vector(*map));
    x = Teuchos::rcp(new Epetra_Vector(*map));

    problem = Teuchos::rcp(new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A, x, b));

    Belos::SolverFactory<double, Epetra_MultiVector, Epetra_Operator> factory;
    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();

    solverParams->set("Maximum Iterations", 1000);
    solverParams->set("Convergence Tolerance", 1.0e-8);

    solver = factory.create("CG", solverParams);

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
    std::vector<double> Values(4);
    std::vector<int> Indices(4);

    int numEntries = 0;
    int local_nx = local_xmax - local_xmin + 1 + 4;

    int i = 0;

    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            Values.clear();
            Indices.clear();

            double c2 = Kx[ARRAY2D(j,k,local_xmin-2,local_ymin-2,local_nx + 1)];
            double c3 = Kx[ARRAY2D(j+1,k,local_xmin-2,local_ymin-2,local_nx + 1)];
            double c4 = Ky[ARRAY2D(j,k,local_xmin-2,local_ymin-2,local_nx + 1)];
            double c5 = Ky[ARRAY2D(j,k+1,local_xmin-2,local_ymin-2,local_nx + 1)];

            numEntries = 0;

            double diagonal = (1+(2*rx)+(2*ry));

            if(1 != k) {
                numEntries++;
                Values.push_back(-ry*c4);
                Indices.push_back(myGlobalIndices_[i]-nx);
            } 
            
            if(ny != k) {
                numEntries++;
                Values.push_back(-ry*c5);
                Indices.push_back(myGlobalIndices_[i]+nx);
            }

            if(1 != j) {
                numEntries++;
                Values.push_back(-rx*c2);
                Indices.push_back(myGlobalIndices_[i]-1);
            }

            if(nx != j) {
                numEntries++;
                Values.push_back(-rx*c3);
                Indices.push_back(myGlobalIndices_[i]+1);
            }

            A->InsertGlobalValues(myGlobalIndices_[i], numEntries, &Values[0], &Indices[0]);
            A->InsertGlobalValues(myGlobalIndices_[i], 1, &diagonal, &myGlobalIndices_[i]);

            i++;
        }
    }

    A->FillComplete();
    //A->Print(std::cout);

    Indices.clear();
    Values.clear();

    i = 0;

    b->PutScalar(0.0);

    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            Indices.push_back(myGlobalIndices_[i]);

            double c2 = Kx[ARRAY2D(j,k,local_xmin-2,local_ymin-2,local_nx + 1)];
            double c3 = Kx[ARRAY2D(j+1,k,local_xmin-2,local_ymin-2,local_nx + 1)];
            double c4 = Ky[ARRAY2D(j,k,local_xmin-2,local_ymin-2,local_nx + 1)];
            double c5 = Ky[ARRAY2D(j,k+1,local_xmin-2,local_ymin-2,local_nx + 1)];

            double value = u0[ARRAY2D(j,k,local_xmin-2, local_ymin-2, local_nx)];

            if(global_ymin == k) {
                value += ry*c4*u0[ARRAY2D(j,k-1,local_xmin-2,local_ymin-2,local_nx)];
            } else if(global_ymax == k) {
                value += ry*c5*u0[ARRAY2D(j,k+1,local_xmin-2,local_ymin-2,local_nx)];
            }

            if(global_xmin == j) {
                value += rx*c2*u0[ARRAY2D(j-1,k,local_xmin-2,local_ymin-2,local_nx)];
            } else if(global_xmax == j) {
                value += rx*c3*u0[ARRAY2D(j+1,k,local_xmin-2,local_ymin-2,local_nx)];
            }

            Values.push_back(value);

            i++;
        }
    }

    b->ReplaceGlobalValues(numLocalElements_, &Values[0], &Indices[0]);
    //b->Print(std::cout);

    Indices.clear();
    Values.clear();

    i = 0;

    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            Indices.push_back(myGlobalIndices_[i]);

            double value = u0[ARRAY2D(j,k, local_xmin-2, local_ymin-2, local_nx)];
            Values.push_back(value);

            i++;
        }
    }

    x->ReplaceGlobalValues(numLocalElements_, &Values[0], &Indices[0]);
    //x->Print(std::cout);

    problem->setOperator(A);
    problem->setLHS(x);
    problem->setRHS(b);

    problem->setProblem();

    solver->setProblem(problem);

    Belos::ReturnType result = solver->solve();

    const int numIters = solver->getNumIters();
    std::cout << "[STEM]: num_iters = " << numIters << std::endl;

    double* updatedValues = new double[numLocalElements_];
    x->ExtractCopy(updatedValues);
    //x->Print(std::cout);

    i = 0;
    for(int k = local_ymin; k <= local_ymax; k++) {
        for(int j = local_xmin; j <= local_xmax; j++) {
            u0[ARRAY2D(j,k,local_xmin-2, local_ymin-2, local_nx)] = updatedValues[i];
            i++;
        }
    }
}
