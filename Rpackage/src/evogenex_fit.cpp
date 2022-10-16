// [[Rcpp::depends(RcppEigen)]]
//
#include "evogenex_hansen.hpp"
#include "evogenex_brown.hpp"

#include <nloptrAPI.h>

// NOTE: There are some function body definition in noptrAPI.h
// Henece it can not be included in different .cpp file as it
// will give multiple definition error for those functions
// Hence we keep the definition of both fitting methods in this file



double myfunc(unsigned n, const double *x, double *grad, void *data)
{
    MLE * mle = static_cast<MLE *>(data);
    return mle->computeLogLik();
}

double mybrown(unsigned n, const double *x, double *grad, void *data)
{
    MLE_Brown * mle = static_cast<MLE_Brown *>(data);
    return mle->computeLogLik();
}


// [[Rcpp::export]]
List brown_fit(int nterm, const IntegerVector nrep,
        const NumericVector &dat,
        const NumericMatrix &bt,
        double gamma)
{
#if KEEP_LOG
    printf("nterm=%d\n", nterm);
    cout << "nrep: " << endl;
    for (int i=0; i<nrep.length(); i++) { 
            cout << " " << nrep[i]; 
    } 
    cout << endl;
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
            Rcpp::Named("gamma") = mle.par[0],
            Rcpp::Named("sigma.sq") = mle.sigmasq,
            Rcpp::Named("loglik") = -0.5*mle.loglik,
            Rcpp::Named("theta") = mle.theta,
            Rcpp::Named("status") = status);
}


// [[Rcpp::export]]
List evogenex_fit(int nterm, 
        const IntegerVector nrep, 
        int nreg,
        const NumericVector &dat,
        const NumericVector &nbranch,
        const NumericVector &beta,
        const NumericVector &epochs,
        const NumericMatrix &bt,
        double alpha, double gamma)
{
#if KEEP_LOG
    printf("nterm=%d\n", nterm);
    cout << "nrep: " << endl;
    for (int i=0; i<nrep.length(); i++) { 
            cout << " " << nrep[i]; 
    }     
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
            Rcpp::Named("loglik") = -0.5*mle.loglik,
            Rcpp::Named("theta") = mle.theta,
            Rcpp::Named("status") = status);
}
