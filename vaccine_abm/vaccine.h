#ifndef VACCINATION_H
#define VACCINATION_H
#include <vector>
#include <cmath>
#include <iostream>

enum class vaccine_type{
    none, 
    pfizer, 
    astrazeneca, 
    moderna}; // Declare an enum for the different tpes of vaccines.

std::ostream &operator<< (std::ostream& strm, const vaccine_type& vac_name);
std::vector<std::vector<double>> setup_vaccine_schedule(std::string schedule_file);

// This is for the vaccine parameters. Pfizer etc.
class vaccine_parameters{
public:
    vaccine_parameters(int dose, vaccine_type , std::vector<double> susceptibility, std::vector<double> transmissibility, std::vector<double> asymptomatic); // This is the constructor of the vaccination.
    
    vaccine_type    get_type();
    int             get_dose();
    std::vector<double>          get_transmissibility();
    double          get_proportion_asymptomatic(const int & age_bracket);
    double          get_susceptibility(const int & age_bracket);
    
private:
    vaccine_type  type;
    int dose_number;
    std::vector<double> susceptibility; // Stratified by dose and age - [0][age] first dose, [1][age] second dose.
    std::vector<double> transmissibility; // Stratified by dose and symptom status. - Dose first symptom status second.
    std::vector<double> proportion_asymptomatic; // Proportion of individuals that are asymptomatic - Stratified by dose and symptoms
    
};

// This is for the individual.
class vaccine{
    public:
//    vaccine(vaccine_type );
    vaccine(vaccine_parameters &, const int & age_bracket);

    // Functions to get vaccine information. 
        vaccine_type get_type();
        double  get_efficacy();
        int     get_dose();
        double  get_transmissibility(double & t, int symptom_status);
        
        double  get_susceptibility(double & t);
        double  get_probability_asymptomatic(double &t);
        double  get_time_of_vaccination();
        void    vaccinate_individual(vaccine_parameters &, const int &, double & t);
        double  get_first_time();
//        void    vaccinate(vaccine_parameters & vax, double & t ); // Vaccinate an individual. Original type might need to match?
    
        // Function that returns the type of vaccine? 
    private:
        // Vaccine parameters - My estimate of what is required.
        vaccine_type type;
        double  probability_asymptomatic;       //  What is the probability that this individual will be asymptomatic at time of exposure. Is a function of vaccination.
        double  old_asymptomatic;               //  This is required for determining the lag.
        double  time_of_vaccination;            // Time of vaccination.
        double  time_of_first_dose = std::nan("1");            // Time of first dose.
        double  susceptibility;                 // Should we have a function for this?
        double  old_susceptibility;             // Previous level of susceptibility, used to determine what happens for an individual in the 14 days before relevant antibody development.
//        double  transmissibility;
//        double  old_transmissibility;
        std::vector<double>  transmissibility;
        std::vector<double>  old_transmissibility;
        int     number_doses;
};


#endif
