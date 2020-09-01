#include <RcppEigen.h>
#include <nloptrAPI.h>
#include <iostream>
#include <limits>
#include <Eigen/Dense>

// NOTE: There are some function body definition in noptrAPI.h
// Henece it can not be included in different .cpp file as it
// will give multiple definition error for those functions
// Hence we keep the definition of both fitting methods in this file

const double X_TOL = sqrt(std::numeric_limits<double>::epsilon());

#define KEEP_LOG 0

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;

class MLE {
public:
    int nterm;
    int nrep;
    int nreg;
    int n;
    double alpha;
    double gamma;
    double sigmasq;
    double logLik;
    std::vector<double> par;
    const NumericVector &nbranch;
    const NumericVector &beta;
    const NumericVector &epochs;
    const NumericVector &bt;
    const Map<VectorXd> &dat;
#if KEEP_LOG
    ofstream out;
    int fcount;
#endif
    MatrixXd W;
    MatrixXd V;
    VectorXd theta;
    MLE(int _nterm,
        int _nrep,
        int _nreg,
        double _alpha,
        double _gamma,
        const NumericVector &_nbranch,
        const NumericVector &_beta,
        const NumericVector &_epochs,
        const NumericMatrix &_bt,
        const NumericVector &_dat
    ): nterm(_nterm), nrep(_nrep), nreg(_nreg),
    n(_nterm * _nrep), alpha(_alpha), gamma(_gamma),
    par{alpha, gamma},
    nbranch(_nbranch),
    beta(_beta),
    epochs(_epochs),
    bt(_bt),
    dat(as<Map<VectorXd> >(_dat)),
#if KEEP_LOG
    out("fast.log"),
    fcount(0),
#endif
    W(n, nreg),
    V(n, n),
    theta(nreg)
    {
    }
    void computeWeights() {
        int index = 0;
        double y[nterm];
        for (int term=0; term<nterm; term++) {
            //cout << "Processing for terminal " << term << endl;
            int nb = nbranch[term];
            //cout << "Num branch " << nb << endl;
            for (int br=0; br<nb; br++) {
                double t = epochs[index] - epochs[index+br];
                y[br] = exp(-alpha*t);
            }
            //cout << "Computing y for terminal " << term << endl;
            for (int br=0; br<nb-1; br++) {
                y[br] = y[br] - y[br+1];
            }
            //cout << "Creating W for terminal " << term << endl;
            for (int reg=0; reg<nreg; reg++) {
                double dotp = 0.0;
                //cout << "Computing dotp  regime " << reg << endl;
                //cout << "beta for regime " << reg << endl;
                for (int br=0; br<nb; br++) {
                    dotp += y[br]*(beta[index+br] == reg ? 1 : 0);
                    //cout << " " << (beta[index+br] == reg ? 1 : 0);
                }
                //cout << endl;
                //cout << "Dotp   " << dotp << endl;
                for (int rep=0; rep<nrep; rep++) {
                    int row = rep + term*nrep;
                    //cout << "row=" << row << endl;
                    W(row,reg) = dotp;
                }
            }
            index += nb;
        }
    }

    void computeCovars() {
        for (int termi=0; termi<nterm; termi++) {
            for (int repk=0; repk<nrep; repk++) {
                for (int termj=0; termj<nterm; termj++) {
                    for (int repl=0; repl<nrep; repl++) {
                        int p = repk + termi*nrep;
                        int q = repl + termj*nrep;
                        double ii = bt[termi + termi*nterm];
                        double jj = bt[termj + termj*nterm];
                        double ij = bt[termi + termj*nterm];
                        double part1 = exp(alpha*(-ii - jj + 2*ij));
                        double part2 = (1-exp(-2*alpha*ij));
                        V(p,q) = part1*part2/(2*alpha);
                        if ((termi == termj) && (repk == repl)) {
                            V(p,q) += gamma;
                        }
                    }
                }
            }
        }
    }

    double computeLogLik()
    {
        alpha = par[0];
        gamma = par[1];
 
#if KEEP_LOG
        fcount++;
        out << "###################### IN ComputeLogLik ########################" << endl;
        out << "n: " << n << endl;
        out << "nrep: " << nrep << endl;
        out << "alpha: " << alpha << endl;
        out << "gamma: " << gamma << endl;
#endif

        computeWeights();
#if KEEP_LOG
        out << "W:\n" << W << endl;
#endif

        computeCovars();
#if KEEP_LOG
        out << "V:\n" << V << endl;
#endif

        auto solver2 = V.colPivHouseholderQr();

        MatrixXd L = V.llt().matrixL();
#if KEEP_LOG
        out << "L:\n" << L << endl;
#endif

        auto solver = L.colPivHouseholderQr();
        auto X = solver.solve(W);
#if KEEP_LOG
        out << "X:\n" << X << endl;
#endif

        auto y = solver.solve(dat);
#if KEEP_LOG
        out << "y:\n" << y << endl;
#endif

        theta = X.bdcSvd(ComputeThinU | ComputeThinV).setThreshold(X_TOL).solve(y);
#if KEEP_LOG
        out << "theta:\n" << theta << endl;
#endif

        auto e = (W*theta - dat);
#if KEEP_LOG
        out << "e:\n" << e << endl;
#endif

        auto zz = solver2.solve(e);
#if KEEP_LOG
        out << "zz:\n" << zz << endl;
#endif

        auto q = e.transpose()*zz;
#if KEEP_LOG
        out << "q:\n" << q << endl;
#endif

        sigmasq = q(0,0)/n;
#if KEEP_LOG
        out << "sigmasq: " << sigmasq << endl;
#endif

        double detv = 2 * solver.logAbsDeterminant(); //toDenseMatrix().diagonal().array().log().sum();
        
        logLik = n*log(2*M_PI) + n*(1+log(sigmasq)) + detv;
#if KEEP_LOG
        out << "logLik: " << logLik << endl;
#endif
        return(logLik);
    }
};

double myfunc(unsigned n, const double *x, double *grad, void *data)
{
    MLE * mle = static_cast<MLE *>(data);
    return mle->computeLogLik();
}


class MLE_Brown {
public:
    int nterm;
    int nrep;
    int n;
    double gamma;
    double sigmasq;
    double logLik;
    std::vector<double> par;
    const NumericMatrix &bt;
    const Map<VectorXd> &dat;
#if KEEP_LOG
    ofstream out;
    int fcount;
#endif
    MatrixXd W;
    MatrixXd V;
    VectorXd theta;
    MLE_Brown(int _nterm,
        int _nrep,
        double _gamma,
        const NumericMatrix &_bt,
        const NumericVector &_dat
    ): nterm(_nterm), nrep(_nrep),
    n(_nterm * _nrep),
    gamma(_gamma),
    par{gamma, 0},
    bt(_bt),
    dat(as<Map<VectorXd> >(_dat)),
#if KEEP_LOG
    out("fast_brown.log"),
    fcount(0),
#endif
    W(MatrixXd::Ones(n, 1)),
    V(n, n),
    theta(1)
    {
    }

    void computeCovars() {
        for (int termi=0; termi<nterm; termi++) {
            for (int repk=0; repk<nrep; repk++) {
                for (int termj=0; termj<nterm; termj++) {
                    for (int repl=0; repl<nrep; repl++) {
                        int p = repk + termi*nrep;
                        int q = repl + termj*nrep;
                        V(p,q) = bt[termi + termj*nterm];
                        if ((termi == termj) && (repk == repl)) {
                            V(p,q) += gamma;
                        }
                    }
                }
            }
        }
    }

    double computeLogLik()
    {
        gamma = par[0];
 
#if KEEP_LOG
        fcount++;
        out << "###################### IN ComputeLogLik ########################" << endl;
        out << "n: " << n << endl;
        out << "nrep: " << nrep << endl;
        out << "alpha: " << alpha << endl;
        out << "gamma: " << gamma << endl;
#endif

        //computeWeights();
#if KEEP_LOG
        out << "W:\n" << W << endl;
#endif

        computeCovars();
#if KEEP_LOG
        out << "V:\n" << V << endl;
#endif

        auto solver2 = V.colPivHouseholderQr();

        // V = L * L^T
        MatrixXd L = V.llt().matrixL();
#if KEEP_LOG
        out << "L:\n" << L << endl;
#endif

        auto solver = L.colPivHouseholderQr();
        // X = L^{-1} * W
        auto X = solver.solve(W);
#if KEEP_LOG
        out << "X:\n" << X << endl;
#endif

        // y = L^{-1} * dat
        auto y = solver.solve(dat);
#if KEEP_LOG
        out << "y:\n" << y << endl;
#endif

        // theta = Least square solution of: X * theta = y
        // that is multiply both sides by X^T: X^T * X * theta = X^T * y
        // that is theta = (X^T * X)^{-1} * y
        theta = X.bdcSvd(ComputeThinU | ComputeThinV).setThreshold(X_TOL).solve(y);
#if KEEP_LOG
        out << "theta:\n" << theta << endl;
#endif

        auto e = (W*theta - dat);
#if KEEP_LOG
        out << "e:\n" << e << endl;
#endif

        auto zz = solver2.solve(e);
#if KEEP_LOG
        out << "zz:\n" << zz << endl;
#endif

        auto q = e.transpose()*zz;
#if KEEP_LOG
        out << "q:\n" << q << endl;
#endif

        sigmasq = q(0,0)/n;
#if KEEP_LOG
        out << "sigmasq: " << sigmasq << endl;
#endif

        double detv = 2 * solver.logAbsDeterminant(); //toDenseMatrix().diagonal().array().log().sum();
        
        logLik = n*log(2*M_PI) + n*(1+log(sigmasq)) + detv;
#if KEEP_LOG
        out << "logLik: " << logLik << endl;
#endif
        return(logLik);
    }
};

double mybrown(unsigned n, const double *x, double *grad, void *data)
{
    MLE_Brown * mle = static_cast<MLE_Brown *>(data);
    return mle->computeLogLik();
}


// [[Rcpp::export]]
List brown_fit(int nterm, int nrep,
        const NumericVector &dat,
        const NumericMatrix &bt,
        double gamma)
{
#if KEEP_LOG
    printf("nterm=%d\n", nterm);
    printf("nrep=%d\n", nrep);
    cout << "gamma: " << gamma << endl; 
    cout << "bt: ";
    for (int i=0; i<nterm; i++) {
        for (int j=0; j<nterm; j++) { cout << " " << bt(i,j); } cout << endl;
    }
#endif

    MLE_Brown mle(nterm, nrep, gamma, bt, dat);

    double lb[2] = { 1e-10, 1e-10 }; 		// lower bounds
    double ub[2] = { 1e+10, 1e+10 }; 		// upper bounds

    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_SBPLX, 1); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_xtol_rel(opt, X_TOL);
    nlopt_set_maxeval(opt, 10000);
    nlopt_set_min_objective(opt, mybrown, &mle);

    double minf; // minimum objective value, upon return
    int status = 0;
    if (nlopt_optimize(opt, &(mle.par[0]), &minf) < 0) {
#if KEEP_LOG
        Rcpp::Rcout << "nlopt failed!" << std::endl;
#endif
    } else {
        status = 1;
#if KEEP_LOG
        {
            Rcpp::Rcout << std::setprecision(5)
                << "Found minimum at f(" << mle.par[0] << "," << mle.par[1] << ") "
                << "= " << std::setprecision(8) << minf
                << " after " << mle.fcount << " function evaluations."
                << std::endl;
        }
#endif
    }
    nlopt_destroy(opt);
    mle.computeLogLik();
    return Rcpp::List::create(
            Rcpp::Named("alpha") = mle.par[1],
            Rcpp::Named("gamma") = mle.par[0],
            Rcpp::Named("sigma.sq") = mle.sigmasq,
            Rcpp::Named("logLik") = -0.5*mle.logLik,
            Rcpp::Named("theta") = mle.theta,
            Rcpp::Named("status") = status);
}


// [[Rcpp::export]]
List evogenex_fit(int nterm, int nrep, int nreg,
        const NumericVector &dat,
        const NumericVector &nbranch,
        const NumericVector &beta,
        const NumericVector &epochs,
        const NumericMatrix &bt,
        double alpha, double gamma)
{
#if KEEP_LOG
    printf("nterm=%d\n", nterm);
    printf("nrep=%d\n", nrep);
    printf("nreg=%d\n", nreg);
    cout << "nbranch: ";
    for (auto v:nbranch) { cout << " " << v; } cout << endl;
    cout << "beta: "; for (auto v:beta) { cout << " " << v; } cout << endl;
    cout << "epochs: "; for (auto v:epochs) { cout << " " << v; } cout << endl;
    cout << "alpha: " << alpha << endl; 
    cout << "gamma: " << gamma << endl; 
    cout << "bt: ";
    for (int i=0; i<nterm; i++) {
        for (int j=0; j<nterm; j++) { cout << " " << bt(i,j); } cout << endl;
    }
#endif

    MLE mle(nterm, nrep, nreg, alpha, gamma, nbranch, beta, epochs, bt, dat);

    double lb[2] = { 1e-10, 1e-10 }; 		// lower bounds
    double ub[2] = { 1e+10, 1e+10 }; 		// upper bounds

    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_SBPLX, 2); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_xtol_rel(opt, X_TOL);
    nlopt_set_maxeval(opt, 10000);
    nlopt_set_min_objective(opt, myfunc, &mle);

    double minf; // minimum objective value, upon return
    int status = 0;
    if (nlopt_optimize(opt, &(mle.par[0]), &minf) < 0) {
#if KEEP_LOG
        Rcpp::Rcout << "nlopt failed!" << std::endl;
#endif
    } else {
        status = 1;
#if KEEP_LOG
        {
            Rcpp::Rcout << std::setprecision(5)
                << "Found minimum at f(" << mle.par[0] << "," << mle.par[1] << ") "
                << "= " << std::setprecision(8) << minf
                << " after " << mle.fcount << " function evaluations."
                << std::endl;
        }
#endif
    }
    nlopt_destroy(opt);
    mle.computeLogLik();
    return Rcpp::List::create(
            Rcpp::Named("alpha") = mle.par[0],
            Rcpp::Named("gamma") = mle.par[1],
            Rcpp::Named("sigma.sq") = mle.sigmasq,
            Rcpp::Named("logLik") = -0.5*mle.logLik,
            Rcpp::Named("theta") = mle.theta,
            Rcpp::Named("status") = status);
}
