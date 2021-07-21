#include "nbinrnd.h"

int nbinrnd::operator()(std::default_random_engine & generator){
    
    double lambda = sample_gamma(generator);
    std::poisson_distribution<int> sample(lambda);

    return sample(generator);
}

nbinrnd::nbinrnd(double r_in,double p_in):r(r_in),p(p_in){
    sample_gamma = std::gamma_distribution<double>(r,p/(1.0-p));
}

double nbinrnd::r_return(){
    return r;
}

double nbinrnd::p_return(){
    return p;
}
