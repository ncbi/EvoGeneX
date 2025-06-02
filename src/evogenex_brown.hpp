#ifndef __EVOGENEX_BROWN_HPP__
#define __EVOGENEX_BROWN_HPP__

#include "evogenex_common.h"

class MLE_Brown {
public:
    int nterm;
    const IntegerVector nrep;
    int n;
    double gamma;
    double sigmasq;
    double loglik;
    std::vector<double> par;
    const NumericMatrix &bt;
    const Map<VectorXd> dat;
    std::vector<int> nrep_sums;
#if KEEP_LOG
    ofstream out;
    int fcount;
#endif
    MatrixXd W;
    MatrixXd V;
    VectorXd theta;
    MLE_Brown(int _nterm,
        const IntegerVector _nrep,
        double _gamma,
        const NumericMatrix &_bt,
        const NumericVector &_dat
    ): nterm(_nterm), nrep(_nrep),
    n(std::accumulate(nrep.begin(), nrep.end(), 0)),
    gamma(_gamma),
    par{gamma, 0},
    bt(_bt),
    dat(as<Map<VectorXd> >(_dat)),
    nrep_sums(nterm),
#if KEEP_LOG
    out("fast_brown.log"),
    fcount(0),
#endif
    W(MatrixXd::Ones(n, 1)),
    V(n, n),
    theta(1)
    {
        std::partial_sum(nrep.begin(), nrep.end(), nrep_sums.begin());
        nrep_sums.insert(nrep_sums.begin(), 0);
    }

    void computeCovars() {
        for (int termi=0; termi<nterm; termi++) {
            for (int repk=0; repk<nrep[termi]; repk++) {
                for (int termj=0; termj<nterm; termj++) {
                    for (int repl=0; repl<nrep[termj]; repl++) {
                        int p = repk + nrep_sums[termi];
                        int q = repl + nrep_sums[termj];
                        V(p,q) = bt(termi, termj);
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
        out << "nrep sums: " << endl;
        for (std::vector<int>::const_iterator it=nrep_sums.begin(); it != nrep_sums.end(); it++) { 
                out << " " << *it; 
        } 
        out << endl;
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
        
        loglik = n*log(2*M_PI) + n*(1+log(sigmasq)) + detv;
#if KEEP_LOG
        out << "loglik: " << loglik << endl;
#endif
        return(loglik);
    }
};

#endif // __EVOGENEX_BROWN_HPP__
