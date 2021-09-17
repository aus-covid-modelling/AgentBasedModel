#ifndef IBM_SIMULATION_H
#define IBM_SIMULATION_H
#include <vector>
#include <random>
#define NEW_COVID_MODEL

//      Forward declare required values. 
class individual; 
class community;
class household;
class facility;
class workplace;

class disease_model{
        public:
        disease_model(std::vector<double> beta_H, std::vector<double> beta_C_in, std::vector<std::vector<double>> contact_matrix_in,std::vector<double> b,std::vector<double> w);
        disease_model(double beta_H,double beta_C, std::vector<std::vector<double>> contact_matrix_in,std::vector<double> b,std::vector<double> w);
        std::vector<double> beta;               //  Household Transmission probability. (Maybe define R0 instead)
        std::vector<double> beta_C;             //  Community transmission probability
        std::vector<std::vector<double>> contact_matrix;        // Contacts per day in each age bracket.
        
        // Distributions! The disease model should initialise some distributions that wil be used!
        // For now lets just use exponential distributions.
        std::exponential_distribution<double>   gen_tau_E;        // Time moving from exposed to infected.
        std::exponential_distribution<double>   gen_tau_S;        // Time from infected to symptoms
        std::exponential_distribution<double>   gen_tau_R;              // Time from infected to recovered.
        std::piecewise_constant_distribution<double>      gen_tau_isolation; // Time before symptoms the individual is isolated.

        //      Define the simulation types here. 
        void covid_model(individual &, std::vector<household> &, community &, const double);            // Covid model with one community- Input an individual and update status.
        void individual_covid(individual &);
        
        // COVID-19 Contact model equations.
        // Age stratified contact model (ASCM)
#ifdef NEW_COVID_MODEL
        double  covid_ascm(std::vector<individual>& , std::vector<household> &, std::vector<std::vector<int>> &, double t0, double t1, double dt, std::vector<size_t> & E, std::vector<size_t> & I, std::vector<size_t> & newly_symptomatic);
        double  covid_one_step_ascm(std::vector<individual>& , std::vector<household> &, std::vector<std::vector<int>> &, double t0, double dt, std::vector<size_t> & E, std::vector<size_t> & I, std::vector<size_t> & newly_symptomatic);
        void    infection_ascm(double t, individual &infected_individual, household & house, std::vector<individual> & residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double dt, std::vector<size_t> & newly_exposed);
        bool    distribution_exposed_update(individual & person, size_t & ind_number, std::vector<size_t> & , double & t);
        bool    distribution_infected_update(individual &person, size_t & ind_number, std::vector<size_t> & , double & t);
#endif
    
        double covid_ascm_R0(std::vector<individual> & residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double t0, double t1, double dt, std::vector<size_t> & E, std::vector<size_t> & I,std::vector<size_t> & newly_symptomatic);
        double  covid_one_step_ascm_R0(std::vector<individual>& residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double t0, double dt, std::vector<size_t> & E, std::vector<size_t> & I,std::vector<size_t> & newly_symptomatic);  

        // Facility - used to implement quarantine.
        void   covid_facility(std::vector<individual>& , facility &,  double t0, double t1, double dt);
        double covid_facility_one_step(std::vector<individual>& , facility &, double t,double dt);
      
        // Seed infection
        void    seed_infection(individual &resident, household &house);
        void    seed_infection(individual &resident);
        void    seed_exposure(individual & resident);
        void    seed_exposure(individual & resident, double & t);

        // Expose infect and recover an individual. 
        inline void expose_individual(individual &);
        inline void infect_individual(individual &);
        inline void recover_individual(individual &);
        inline void susceptible_individual(individual & resident);
    
        // Expose infect and recover an individual.
        inline void expose_individual(individual &, double &);
        inline void infect_individual(individual &, double &);
        inline void recover_individual(individual &, double &);
        inline void susceptible_individual(individual & resident, double &);
};
#endif
