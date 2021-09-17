#ifndef NBINRND_H
#define NBINRND_H
#include<random>
class nbinrnd
{
public:

    // constructor and reset functions
    nbinrnd(double k = 1, double p = 0.5);

    // Parameters
    double r; 
    double p;
    std::gamma_distribution<double> sample_gamma;

    int operator()(std::default_random_engine &);

    double r_return();
    double p_return();

};
#endif
