#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include "../vaccine_abm/individual.h"
#include "../vaccine_abm/households.h"
#include "../vaccine_abm/ibm_simulation.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <algorithm>


// Define the information required for the contact ibm.
static std::random_device rd;
static std::default_random_engine generator(rd());
static std::uniform_real_distribution<double> genunf_std(0.0,1.0);


// Define the structure of parameters.
struct parameter_struct{
  //  Disease model parameters. (Will be varied)
    double  beta; // Household beta.
    int     num_houses;
    double  average_house_size;
};

struct events{
public:
    double threshold_50_time = -10.0;
    double threshold_60_time = -10.0;
    double threshold_70_time = -10.0;
    double threshold_80_time = -10.0;
};

std::vector<individual> run_model(double beta_C, parameter_struct parameters, std::vector<double> & age_brackets,std::vector<double> & vaccinated_proportion, std::vector<double> & population_pi, std::vector<std::vector<double>> &contact_matrix, vaccine_parameters & no_vaccine, std::vector<vaccine_parameters> & pfizer, std::vector<vaccine_parameters> & astrazeneca, std::vector<vaccine_parameters> & moderna,std::vector<std::vector<double>> & tti_distribution){
    // Inputs from QUANTIUM - pfizer doses per week etc are the proportion of doses per week in each age bracket.
    
//    std::vector<double> b{60,61};
//    std::vector<double> w{1};
    // Number of age brackets.
    size_t num_brackets = population_pi.size();

    //  Disease model parameters. (Will be varied)
    double beta  = parameters.beta; // Household infection rate.
    
    // Disease model initialisation.
    disease_model covid(beta, beta_C, contact_matrix,tti_distribution[0],tti_distribution[1]); // Load the disease model, asymptomatic infections and severity are associated with each individual.
    
    
    //  Parameters for households.
    int     num_houses              = parameters.num_houses;
    double  average_house_size      = parameters.average_house_size;

    // Storage for households and residents.
    std::vector<household> houses; houses.reserve(num_houses);
    std::vector<individual> residents; residents.reserve(std::floor(num_houses*average_house_size));

    // Initialise all of the community...
    community city(0,num_houses,average_house_size,residents,houses,age_brackets,population_pi,no_vaccine); // Vaccinated proportion not used here now.
    std::cout << "residents = " << residents.size() << " houses = " << houses.size() << "\n";

    // Vaccine Schedule - stratified by age group, randomly permute individuals into the order that they will get the first dose.
    std::vector<std::vector<int>> age_dependent_vaccination(9,std::vector<int>(0,0)); // List the individuals that will be vaccinated.
    for(int i = 0; i < residents.size(); i++){                  // For each individual.
        if(residents[i].age<16.0){ continue;}                   // You are not elligible.
            int row_ref = std::floor(residents[i].age/10.0);    // Age bracket.
            row_ref = (row_ref>8)?8:row_ref;                    // Which vaccination strata are you in.
            age_dependent_vaccination[row_ref].push_back(i);
    }

    // Timestepping parameters
    double t = 0.0;
    double t_end = 60; // 6 months past the end of the simulation.
    double dt = pow(2,-4.0); // Probably doesnt need 16 per day but it runs so quick...
    
    // Create memory that tracks who is exposed, E_ref, and who is infected, I_ref. gen_res is used to sample from the list of residents uniformly.
    std::uniform_int_distribution<int> gen_res(0,(int) residents.size()-1); // Dont use rand() its apparently bad (something something entropy)
    std::vector<size_t> E_ref; E_ref.reserve(10000); // Magic number reserving memory.
    std::vector<size_t> I_ref; I_ref.reserve(10000); // Magic number of reserved.
    
    // Assemble the age_matrix (this is a list of people that are in each age bracket).
    std::vector<std::vector<int>> age_matrix(num_brackets);
    assemble_age_matrix(residents,age_matrix); // Nobody moves from the age matrix so only have to do it once.
    
    int cluster_ref = 0; // Track the exposure number, can show that one dominates.
        while(cluster_ref < 1){
            int exposed_resident = gen_res(generator); // Randomly sample from all the population.
                if(residents[exposed_resident].vaccine_status.get_type()==vaccine_type::none){
                    if(residents[exposed_resident].covid.infection_status!='E'){
                    covid.seed_exposure(residents[exposed_resident],t); // Random resident has become infected
                    residents[exposed_resident].covid.cluster_number = cluster_ref;
                    cluster_ref ++ ; // Increment the number of clusters.
                    E_ref.push_back(exposed_resident); // Start tracking them.
        }
    }
        }
    
    while(t < t_end){
        // Infection model!
        std::cout << "Current time = " << t << "\n";
        // Reinitialise the new symptomatic infections.
        std::vector<size_t> newly_symptomatic; newly_symptomatic.reserve(3000); // Reserve size so that reallocation isnt neccesary. Its a magic number.
        
            t = covid.covid_ascm_R0(residents,houses,age_matrix,t,t+dt,dt,E_ref,I_ref,newly_symptomatic);
    
    }
    // Finished model.
    return residents;
}

// The main script driver that runs and call run_model.
int main(int argc, char *argv[]){

    if(argc == 1){
        std::cout << "ERROR: Parameter file not loaded!" << std::endl;
        return 0;
    }
    // Parameter read in and simulation here!
    std::string filename(argv[1]);
    std::string file_type(".txt");
    
    try{
        if(filename.compare(filename.size()-4,file_type.size(),file_type)!=0){
            throw std::logic_error("Filename must end in .txt");
        }
    }catch(std::logic_error &e){
        filename = filename + file_type;
        std::cout << "Appending .txt to filename... filename is " << filename << "\n";
    }

    // Read in the parameter files.
    std::ifstream parameters_in(filename);
    std::string line;
    std::string var_name;
    std::string tp_filename;
    std::string scenario_ref;
    std::string tti_filename;
    std::string SEED_INFECTION;
    // Create folder
    std::string folder;

    // Disease model parameters structure
    parameter_struct parameters;
    int num_sims = 0;

    // Read in the parameters file.
    while(!parameters_in.eof()){
        parameters_in >> var_name;
        // Assign to appropriate parameters.
        if(var_name.compare("num_sims")==0){
            parameters_in >> num_sims;
            std::cout << num_sims << std::endl;
        }if(var_name.compare("Beta_H")==0){
            parameters_in >> parameters.beta;
            std::cout << parameters.beta << std::endl;
        }else if(var_name.compare("Folder_name")==0){
            parameters_in >> folder;
        }else if(var_name.compare("num_houses")==0){
            parameters_in >> parameters.num_houses;
        }else if(var_name.compare("average_house_size")==0){
            parameters_in >> parameters.average_house_size;
        }else if(var_name.compare("TransmissionPotential")==0){
            parameters_in >> tp_filename;
        }else if(var_name.compare("VaccineScenario")==0){
            parameters_in >> scenario_ref;
        }else if(var_name.compare("TTIDistribution")==0){
            parameters_in >> tti_filename;
        }else if(var_name.compare("SeedInfection")==0){
            parameters_in >> SEED_INFECTION;
        }
    }
    parameters_in.close(); // Close the file.
    
    
    // We have a vector of tranmission potentials we want to sample from.
    std::ifstream TP_read(tp_filename); // Load posterior of transmission potentials.
    std::vector<double> TP(0,0); // Should be read in from Freya and Nic
    if(TP_read.is_open()){
    double temp_TP;
        while(TP_read >> temp_TP){
            TP.push_back(temp_TP);
        }
    TP_read.close();
    }else{
        throw std::logic_error("File TP_posterior.txt does not exist in the working directory.");
    }
    
    if(num_sims > (int) TP.size()){
        num_sims = (int) TP.size();
    } // Ensure that we dont try and simulate from more than the number of posterior samples. Otherwise we would not be sampling right.
    
    // Which values of TP will we sample from.
    std::vector<size_t> TP_ref(TP.size(),0);
    std::iota(TP_ref.begin(),TP_ref.end(),0);
    std::shuffle(TP_ref.begin(),TP_ref.end(),std::mt19937{std::random_device{}()}); // Randomly shuffle the vector. We will be sample from the first num_sims.
    
    // Create folder.
    std::string directory = "./outputs/" + folder;
    int main_folder = mkdir(directory.c_str(),0777); // Create folder.
    (void) main_folder; // Unused variable;
    

    //  Population demographic parameters (There should be 16 here)
    std::vector<double> population_pi{0.061,0.063,0.061,0.057,0.066,0.074,0.075,0.071,0.063,0.065,0.059,0.061,0.056,0.049,0.043,0.031+0.021+0.022};
    int num_brackets    = (int) population_pi.size();
    std::vector<double> vaccinated_proportion(num_brackets,0.0);
    
    if(vaccinated_proportion.size()!=population_pi.size()){
        throw std::logic_error("Vaccinated proportion and population stratification do not match in size.");
    }
    
    // Define the age brackets that are used in the contact_matrix (this is used to generate an age for each individual)
    std::vector<double> age_brackets{0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75};
    if(age_brackets.size()!=num_brackets){
        throw std::logic_error("Error: age brackets not the same size as population proportion!!");
    }
    
    // Contact matrix (currently importing Prem's Australia matrix).
    std::vector<std::vector<double>> contact_matrix(num_brackets,std::vector<double>(num_brackets,0.0));
    std::ifstream matrix_read;
    matrix_read.open("contact_matrix.txt");
    if(matrix_read.is_open()){
    for(size_t i = 0; i < num_brackets;i++){
        for(size_t j = 0; j <num_brackets;j++){
            matrix_read >> contact_matrix[i][j];
        }
    }
    matrix_read.close();
    }else{
        throw std::logic_error("contact_matrix.txt not found in working directory.\n");
    }
   
    
    std::vector<std::vector<double>> tti_distribution;
    std::ifstream tti_stream(tti_filename + ".csv");

    if(tti_stream.is_open()){
        std::string line;
        double value;
        while(std::getline(tti_stream,line)){
            std::stringstream stream_line(line);
            std::string row_val;
            std::vector<double> row;
            while(std::getline(stream_line,row_val,',')){
                std::stringstream stream_row(row_val);
                stream_row >> value;
                row.push_back(value);
            }
            tti_distribution.push_back(row);
        }
        tti_stream.close();
    }else{
        throw std::logic_error("The tti distribution file was not found.\n");
    }
    
    
    // Vaccine parameters. - DELTA
    // No vaccination.
    std::vector<double> tau_none{1,0.5}; // Currently symptomatc/asymptomatic.
    std::vector<double> q_none{0.29, 0.29, 0.21, 0.21, 0.27, 0.27, 0.33, 0.33, 0.4, 0.4, 0.49, 0.49, 0.63, 0.63, 0.69, 0.69}; // This is symptomatic not asymptomatic.
    std::vector<double> xi_none{0.4, 0.4, 0.38, 0.38, 0.79, 0.79, 0.86, 0.86, 0.8, 0.8, 0.82, 0.82, 0.88, 0.88, 0.74, 0.74};
    vaccine_parameters no_vaccine(0,vaccine_type::none, xi_none, tau_none, q_none); // No vaccination will be passed into everyone. Required for the start of the simulation.
    
    // Pfizer vaccination parameters. (I would like to have all of this read in)
    std::vector<double> tau_P_1{0.67,0.335};;
    std::vector<double> q_P_1{0.1943, 0.1943, 0.1407, 0.1407, 0.1809, 0.1809,0.2211, 0.2211, 0.268, 0.268, 0.3283, 0.3283, 0.4221, 0.4221, 0.4623, 0.4623};
    std::vector<double> xi_P_1{0.28,0.28,0.266,0.266,0.553,0.553,0.602,0.602,0.56,0.56,0.574,0.574,0.616,0.616,0.518,0.518};
    std::vector<double> tau_P_2{0.17, 0.085};
    std::vector<double> q_P_2{0.0493,0.0493,0.0357,0.0357,0.0459,0.0459,0.0561,0.0561,0.068,0.068,0.0833,0.0833,0.1071,0.1071,0.1173,0.1173};
    std::vector<double> xi_P_2{0.084,0.084,0.0798,0.0798,0.1659,0.1659,0.1806,0.1806,0.168,0.168,0.1722,0.1722,0.1848,0.1848,0.1554,0.1554};
    std::vector<vaccine_parameters> pfizer{vaccine_parameters(1,vaccine_type::pfizer,xi_P_1,tau_P_1,q_P_1),vaccine_parameters(2,vaccine_type::pfizer,xi_P_2,tau_P_2,q_P_2)};

    // Astrazeneca vaccination parameters.
    std::vector<double> tau_AZ_1{0.67,0.335};
    std::vector<double> q_AZ_1{0.1943,0.1943,0.1407,0.1407,0.1809,0.1809,0.2211,0.2211,0.268,0.268,0.3283,0.3283,0.4221,0.4221,0.4623,0.4623}; // Probability of being symptomatic.
    std::vector<double> xi_AZ_1{0.328,0.328,0.3116,0.3116,0.6478,0.6478,0.7052,0.7052,0.656,0.656,0.6724,0.6724,0.7216,0.7216,0.6068,0.6068}; // Suscseptibility.
    std::vector<double> tau_AZ_2{0.39, 0.195};
    std::vector<double> q_AZ_2{0.1131,0.1131,0.0819,0.0819,0.1053,0.1053,0.1287,0.1287,0.156,0.156,0.1911,0.1911,0.2457,0.2457,0.2691,0.2691};
    std::vector<double> xi_AZ_2{0.16,0.16,0.152,0.152,0.316,0.316,0.344,0.344,0.32,0.32,0.328,0.328,0.352,0.352,0.296,0.296};
    std::vector<vaccine_parameters> astrazeneca{vaccine_parameters(1,vaccine_type::astrazeneca,xi_AZ_1,tau_AZ_1,q_AZ_1),vaccine_parameters(2,vaccine_type::astrazeneca,xi_AZ_2,tau_AZ_2,q_AZ_2)};

    // Moderna vaccination parameters. (I would like to have all of this read in) - assumed to be the same as Pfizer.
    std::vector<vaccine_parameters> moderna{vaccine_parameters(1, vaccine_type::moderna, xi_P_1, tau_P_1, q_P_1),vaccine_parameters(2, vaccine_type::moderna, xi_P_2, tau_P_2, q_P_2)};

    // The R0 calculation that has been updated to account for all of the heterogeneity built in at the moment....
    double tau_s = tau_none[0]; double tau_sc = tau_none[1];
    double sum_expression = 0.0;
    for(int k = 0; k < (int) num_brackets;k++){
        double xi_k = no_vaccine.get_susceptibility(k);
        
        double internal_sum = 0.0;
        for(int i = 0; i < (int) num_brackets; i++){
            double lambda_ik = contact_matrix[i][k];
            internal_sum += lambda_ik*((tau_s - tau_sc)*no_vaccine.get_proportion_asymptomatic(i) + tau_sc)*population_pi[i];
        };
        sum_expression += internal_sum*xi_k;
    };
    
    // Beta_C is calculated in a for loop because we are samping over multiple Transmission potentials.
    std::cout << num_sims << " simulations " << std::endl;
    for(int i = 0; i < num_sims; i++){
        // Calculate beta_C for the simulation.
        double beta_C = TP[TP_ref[i]]/((5.0)*sum_expression);
        std::cout << "Sim Number = " << i+1 << "\n";
        std::cout << "R0 = " << TP[TP_ref[i]] << "\n";
        std::cout << "beta C = " << beta_C << std::endl;
        auto begin = std::chrono::high_resolution_clock::now();
        std::vector<individual> residents = run_model(beta_C, parameters, age_brackets, vaccinated_proportion, population_pi, contact_matrix, no_vaccine, pfizer, astrazeneca, moderna,tti_distribution);
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
        
        // Filename.
        std::string filename = directory  + "/sim_number_" + std::to_string(i+1) + ".csv";
        // Write file.
        std::ofstream output_file(filename);

        if(output_file.is_open()){

            output_file << "Individual, Age, Current Vaccine, Current doses, Time of first dose, Time of last dose, Vaccine at infection, Doses at infection, Time latest dose at infection, Severity, Symptomatic, Time isolated, Detected case, Time of detection, Time of exposure, Time of infection, Time of symptom onset, Time of recovery, Cluster, Secondary Infections" << std::endl;
            
            int ind_num = 0;
            // To do: we must define both the vaccination status of the individual at infection and at end?
            for(individual & person: residents){
                bool is_infected = !std::isnan(person.covid.time_of_exposure); // Are they infected.
                vaccine_type  vac = person.vaccine_status.get_type();
                std::string vaccine_name = (vac==vaccine_type::none)?"None":(vac==vaccine_type::pfizer)?"Pfizer":(vac==vaccine_type::astrazeneca)?"Astrazeneca":"Moderna";
                
                // Infected write!
                if(is_infected){
                    // Determine vaccine at time of infection.
                    vaccine_type infection_vac = person.infection_statistics.vaccine_status;
                    std::string vac_at_infection = (infection_vac==vaccine_type::none)?"None":(infection_vac==vaccine_type::pfizer)?"Pfizer":(infection_vac==vaccine_type::astrazeneca)?"Astrazeneca":(infection_vac==vaccine_type::moderna)?"Moderna":"error";
                    
                    output_file << ind_num << ", " << person.age << ", ";
                    
                    output_file << vaccine_name << ", " << person.vaccine_status.get_dose() <<", " << person.vaccine_status.get_first_time() << ", " << ((person.vaccine_status.get_dose()==2)?person.vaccine_status.get_time_of_vaccination():std::nan("1")) << ", "; // Vaccination at end of simulation.
                    
                    output_file << vac_at_infection << ", " << person.infection_statistics.doses << ", " << person.infection_statistics.time_of_last_dose << ", ";
                    
                    output_file << (person.covid.severe?"Severe":"Mild") << ", " << (person.covid.asymptomatic?"Asymptomatic":"Symptomatic") << ", "<< person.time_isolated << ", " << (person.covid.detected_case?"Detected":"Undetected") << ", " << person.covid.time_of_detection << ", ";
                    
                    output_file << person.covid.time_of_exposure << ", " << person.covid.time_of_infection << ", " << person.covid.time_of_symptom_onset << ", ";
                    output_file << person.covid.time_of_recovery << ", " << person.covid.cluster_number << ", " << person.who_infected.size() << std::endl;
                }else{

                }
                ind_num ++;
            }


        // Close file.
        output_file.close();

        };

        }
    
    return 1;
 
};



