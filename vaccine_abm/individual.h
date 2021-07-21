// Guard for individual class
# ifndef INDIVIDUAL_H
# define INDIVIDUAL_H
#include <vector>
#include <iostream>
#include<random>
#include "vaccine.h" // Require the vaccine header. 

class household;        //  Class definition for households.
class test;
class vaccine; 

// Define a contact class. This is just a std::pair... silly me.
class trace_contact{
    public:
    trace_contact();
    trace_contact(int ind_number, double t);
    void update_time(double);
    double time;
    int who; 
       // Overload the equals operator so as to find the resident number.
        bool operator==(const int & obj2) const{
        if(this->who==obj2){
            return true;
        }else{
            return false;
        }
    };
};

class simulation_statistics
{
public:
    simulation_statistics(); // Trivial constructor.
    void set_infection_stats(vaccine_type , int doses, double time_latest);
    vaccine_type vaccine_status; // This will be for at the time of infection.
    int doses;
    double time_of_last_dose = std::nan("1");
    
private:

    
};

class disease
{
    public:
        //  Constructor functions
        disease(char );
        // Public member data
        char    infection_status;               //  Determine if individual is S E I R or whatever you want really... 
        bool    symptoms = false;               //  Does the individual have symptoms. This will be based off some probability. 
        bool    asymptomatic = false;           //  Is the individual asymptomatic. 
        bool    severe = false;                 //  Is this a severe disease.
        
        // Vaccination parameters... This is a place holder as I debate if this is the location best suited or not.

        // Variables that are assigned at time of exposure.
        double  transmissibility = std::nan("1"); // Assigned at time of exposure. Is a function of vaccination. This is a constant throughout transmission so reassigned here at exposure (else just use vaccine::get_transmissibility(t).
        
        // When does event occur.
        double  time_of_exposure         = std::nan("1");
        double  time_of_symptom_onset    = std::nan("2");
        double  time_of_infection        = std::nan("3");
        double  time_of_recovery         = std::nan("4");
    
        // Latch to determine if pre-symptomatic time is over.
        bool    check_symptoms           = true; // Symptomatic latch.
        
        // Case detections.
        bool    detected_case           = false; // Were they a detected case of COVID.
        double  time_of_detection       = std::nan("5");

        // Future development.
         int cluster_number = -1; // We can track the clusters through time. That could be fun. It can be passed from exposure to exposure.
};

//  Define class for each individual in poopulation.
class individual
{
public:

    //  Constructor functions.
//    individual(double age, int age_bracket_in, int home_id, int community_id, vaccine_type type); // Vaccination constructor
    individual(double age, int age_bracket_in, int home_id, int community_id, vaccine_parameters & Vaccination); // Vaccination constructor for an individual.
    
    //  Parameters of individual class.
    int     age_bracket;                    //  Individual age bracket
    double  age;                            //  Person's age
    int     home_id;                        //  What house do you live in (this is your primary residence). 
    int     community_id;                   //  What community do you belong to. 
    int     current_home;                   //  Where are they at this point in time.
    int     current_community;              //  Where are they for this timestep.

    // Statistics.
    simulation_statistics infection_statistics;
    //  Disease information. (Should we look at multiple strains? What about waning immunity?)
    disease covid;

    // Vaccinated parameters. (Do we want to make a vaccine class.. probably)
    vaccine vaccine_status;
    
    // Intervention parameters - Testing
    bool    tested      = false;            //  Has the individual been tested. 
    double  time_until_test;              // Time until this individual will report for testing?

//    double  time_since_last_test = 100000.0;    // Time since they were last tested
//    bool    awaiting_test_result = false;          // Probably not required for vacinne
//    bool    received_positive_test = false;        // Probably not required for cvaccination.
    
    // Intervention parameters - Isolation and quarantine.
    // Included for inter-provincial travel quarantine measures. 
    bool    quarantined = false;            //  Included for the case when people are quarantined for reasons external to their household. 
    bool    external_quarantine = false;    // Are they externally quarantined

    // Is the individual isolated.
    bool    isolated = false;                   //  Is the individual isolated (because they have become sick)
    bool    external_isolation = false;         // Are they to be isolated away from everyone.
    double  time_isolated = std::nan("4");;             // Time they are in isolation until.
    
    // Could remove these for memory issues.
    //    std::pair<double , int> contacts; // This is the more appropriate way.
    std::vector<trace_contact> contacts; 
    std::vector<trace_contact> household_contacts;  // People that have lived in the house.
    std::vector<int> who_infected;                  //  Who did they infect.

    //  Functions to update individual.
    void update_age(const double & dt);                         //  Update the age of an individual     

};

// Age_sort and overloaded variations.
int age_sort(individual & , std::vector<double>);
int age_sort(double age,std::vector<double> age_brackets);
double generate_age(int bracket, const std::vector<double> & age_brackets); // Generates an age between the chosen age bracket.

# endif
