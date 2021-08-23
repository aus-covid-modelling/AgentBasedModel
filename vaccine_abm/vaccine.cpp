#include "vaccine.h"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
// In this code, it is assumed that it has been long enough to get your previous doses immunity!

//Serialisation for vaccine types
//If new vaccines are added then they will need to go in here.
std::ostream &operator<< (std::ostream& strm, const vaccine_type& vac_name) {
    
    switch (vac_name) {
        case vaccine_type::none:
            strm << "None";
            break;
        case vaccine_type::pfizer:
            strm << "Pfizer";
            break;
        case vaccine_type::moderna:
            strm << "Moderna";
            break;
        case vaccine_type::astrazeneca:
            strm << "AstraZeneca";
            break;
        default:
            strm << "Unknown";
    }
    return strm;
}

vaccine_parameters::vaccine_parameters(int dose, vaccine_type vac_name, std::vector<double> xi_in, std::vector<double> tau_in, std::vector<double> p_a):type(vac_name),susceptibility(xi_in),transmissibility(tau_in),proportion_asymptomatic(p_a),dose_number(dose){
    // Throw error in constructor if the dimensions are off.
    if(p_a.size()!=xi_in.size()){
        throw std::logic_error("xi_in and p_a do not have the same size in the vaccine parameters constructor.");
    }
    
    
}; // This is the constructor of the vaccination.

vaccine_type vaccine_parameters::get_type(){return type;};

int vaccine_parameters::get_dose(){return dose_number;};

/**
 * Set up vaccine schedule from CSV. Convenience function to read in the vaccine schedule
 * from a given CSV.
 * @param schedule_file String indicating path to the schedule file
 * @return Vector containing the number of doses per age-group per week.
 */
std::vector<std::vector<double>> setup_vaccine_schedule(std::string schedule_file) {
    
    std::ifstream schedule_ifstream(schedule_file);
    std::vector<std::vector<double>> schedule_doses_per_week;

    if(schedule_ifstream.is_open()){
        std::string line;
        double value;
        while(std::getline(schedule_ifstream,line)){
            std::stringstream stream_line(line);
            std::string row_val;
            std::vector<double> row;
            while(std::getline(stream_line,row_val,',')){
                std::stringstream stream_row(row_val);
                stream_row >> value;
                row.push_back(value);
            }
            schedule_doses_per_week.push_back(row);
        }
        schedule_ifstream.close();
    } else {
        throw std::logic_error("The schedule file for" + schedule_file + "was not found in vaccination_input.\n");
    }

    return (schedule_doses_per_week);
};

double vaccine_parameters::get_susceptibility(const int & age_bracket){
    if(age_bracket>=(int)susceptibility.size()){
        throw std::logic_error("Error in vaccine_parameters::get_susceptibility(int age_bracket). Age bracket greater than size of susceptibility vector.");
    }
    return susceptibility[age_bracket]; // Dose has a - 1
};

std::vector<double> vaccine_parameters::get_transmissibility(){return transmissibility;};

double vaccine_parameters::get_proportion_asymptomatic(const int & age_bracket){return proportion_asymptomatic[age_bracket];};

// Vaccine constructor for the individual. This is different to the vaccinate function.
vaccine::vaccine(vaccine_parameters & vaccination, const int & age_bracket):type(vaccination.get_type()){
    
    double t = -20.0; // This is the constructor, so make it before time = 0.0; It is assumed that you have this immunity levels. 
    time_of_vaccination = t;
    number_doses = vaccination.get_dose(); // Initialised with the vaccine_parameters information.
   
    // Whatever parameters that are required should be used here. Order does not matter as it is being initialised....which by definition is unvaccinated?
    susceptibility = vaccination.get_susceptibility(age_bracket);
    old_susceptibility = susceptibility;
    
    
    transmissibility = vaccination.get_transmissibility();
    old_transmissibility  = transmissibility;
    
    
    probability_asymptomatic = vaccination.get_proportion_asymptomatic(age_bracket);
    old_asymptomatic = probability_asymptomatic;
    

};

vaccine_type vaccine::get_type(){
    return type;
}

int vaccine::get_dose(){
    return number_doses;
}

double vaccine::get_transmissibility(double & t, int symptom_status){
    // t - this is the current time.
    return ((t - time_of_vaccination)>14.0)*(transmissibility[symptom_status]-old_transmissibility[symptom_status]) + old_transmissibility[symptom_status]; // It is a switch that turns on immune response to vaccination. Should record time of vaccination and infection?
};

double vaccine::get_susceptibility(double & t){
    // t - this is the current time.
    return ((t - time_of_vaccination)>14.0)*(susceptibility-old_susceptibility) + old_susceptibility; // Could have time dependent... Takes 14 days to get immunity effects. Woo lag.
};

double vaccine::get_probability_asymptomatic(double & t){
    // t - this is the current time.
    return ((t - time_of_vaccination)>14.0)*(probability_asymptomatic-old_asymptomatic) + old_asymptomatic; // Could have time dependent... Takes 14 days to get immunity effects. Woo lag.
};

double vaccine::get_time_of_vaccination(){
    return time_of_vaccination;
}

double vaccine::get_first_time(){
    return time_of_first_dose;
}

void vaccine::vaccinate_individual(vaccine_parameters & vaccination, const int & age_bracket, double & t){
    
    if(type!=vaccine_type::none){
        if(type != vaccination.get_type()){
            throw std::logic_error("Vaccination types do not match");
        }
    }
    
    if(number_doses+1!=vaccination.get_dose()){
        throw std::logic_error("Dosage number does not match");
    }
    
    type = vaccination.get_type();
    time_of_vaccination = t;
    number_doses = vaccination.get_dose(); // Initialised with the vaccine_parameters information.
    
    if(number_doses==1){
        time_of_first_dose = t;
    }
    // Whatever parameters that are required should be used here
    
    // Order matters here. old first. 
    old_susceptibility = susceptibility;
    susceptibility = vaccination.get_susceptibility(age_bracket);
    

    old_transmissibility  = transmissibility;
    transmissibility = vaccination.get_transmissibility();
    
    old_asymptomatic = probability_asymptomatic;
    probability_asymptomatic = vaccination.get_proportion_asymptomatic(age_bracket);
    

}
