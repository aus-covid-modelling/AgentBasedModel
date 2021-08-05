#include <chrono>
#include "individual.h"
#include "vaccine.h"
static std::random_device rd;

// Allows for a different random seeding method on Windows. _WIN32 should be defined even on 64 bit.
#ifdef _WIN32
    static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
#else
    static std::default_random_engine generator(rd());
#endif
static std::uniform_real_distribution<double> genunf_std(0.0,1.0);

//  Define constructor for the disease class (removed the trivial constructor)
disease::disease(char status):infection_status(status){}

// Individual constructor - modified for the vaccination_parameters class.
individual::individual(double age, int age_bracket_in, int home_id, int community_id, vaccine_parameters & Vaccination):age_bracket(age_bracket_in),age(age),home_id(home_id),community_id(community_id),covid(disease('S')),vaccine_status(vaccine(Vaccination,age_bracket_in)){
    current_community =     community_id;   // Assign to home community
    current_home      =     home_id;        // Assign to household.
}

void individual::update_age(const double &dt){
    age += dt/365.0;  //  Age is in years and dt is in days.
//    time_since_last_test += dt;
}

double generate_age(int bracket, const std::vector<double> &age_brackets){
   // This function has hardcoded the oldest is the age bracket + 10. 
   // Uniformly samples age between the two brackets. 
    double r = genunf_std(generator);
    double diff;

    if(bracket == (int) age_brackets.size()-1){
        diff = 10; 
    }else{
        diff = age_brackets[bracket+1] - age_brackets[bracket];
    }
    
    return age_brackets[bracket] + r*diff;  
}

int age_sort(individual &person,std::vector<double> age_brackets){
    // Return the age bracket of the individual. Ths function is hardcoded to be in groups of 10 years.
    // Agesort needs a minus one for the final one because referenced from zero!
    double age = person.age;
    for(int i = 0; i < (int) age_brackets.size()-1; i++){
        
        if(age < age_brackets[i+1]){
            person.age_bracket = i;
            return i;
        }
    };  
    person.age_bracket = (int) age_brackets.size()-1;
    return (int) (age_brackets.size()-1);
}

int age_sort(double age,std::vector<double> age_brackets){
    // Return the age bracket of the individual. Ths function is hardcoded to be in groups of 10 years.
    // Agesort needs a minus one for the final one because referenced from zero!
    for(int i = 0; i < (int) age_brackets.size()-1; i++){
        
        if(age < age_brackets[i+1]){
            return i;
        }


    };

    return (int) age_brackets.size()-1;
}

// Contact class. 
trace_contact::trace_contact(int ind_number, double t):time(t),who(ind_number){}
trace_contact::trace_contact(){}

void trace_contact::update_time(double t){
    time = t;
}

simulation_statistics::simulation_statistics(){};
void simulation_statistics::set_infection_stats(vaccine_type vax, int number_doses, double time_latest){
    vaccine_status = vax;
    doses = number_doses;
    time_of_last_dose = time_latest;
}
