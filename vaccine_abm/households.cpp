#include "individual.h"
#include "households.h"
#include <algorithm>
#include <chrono>
#include "vaccine.h"

static std::random_device rd;

// Allows for a different random seeding method on Windows. _WIN32 should be defined even on 64 bit.
#ifdef _WIN32
    static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
#else
    static std::default_random_engine generator(rd());
#endif
static std::uniform_real_distribution<double> genunf_std(0.0,1.0);

// Constructor for community class with age brackets
community::community(int community_id,int num_houses, double house_average, std::vector<individual> & residents,
     std::vector<household> & houses, const std::vector<double> &age_brackets, const std::vector<double> pi, vaccine_parameters & initial_vax_status){
        
        std::poisson_distribution<int> generate_household_size(house_average-1.0); // Subtract 1 and add 1 later? 
        for(int i = 0; i < num_houses; i++){
            //  How many residents in house. 
            int num_residents = 0;
            while(num_residents==0){
            num_residents   = generate_household_size(generator)+1.0;    //  Number of residents in each house. 
            }
            int house_ID        = (int) houses.size();    //  ID for houses to calculate efficient force of infection.
            houses.push_back(household(community_id, house_ID, num_residents, residents, age_brackets,pi,initial_vax_status));   //  Initialise the houses.
            population_size += num_residents;
        }

}

//  Household with age_brackets.
household::household(int community_id, int house_ID, int num_residents, std::vector<individual> & residents, const std::vector<double> & age_brackets,const std::vector<double> pi,vaccine_parameters & initial_vax_status){
    // When a household is constructed, we must construct the vector of residents that live there
    std::discrete_distribution<int> gen_age_bracket(pi.begin(),pi.end());
    for(int i = 0; i < num_residents; i++){
        int bracket = gen_age_bracket(generator);
        double age = generate_age(bracket, age_brackets);
        current_residents.push_back((int) residents.size());
        residents.push_back(individual(age,bracket,house_ID,community_id,initial_vax_status));
    };
}


void household::add_resident(int resident_number){
    // Add resident to the current_residents.
    current_residents.push_back(resident_number);
}

void household::remove_resident(int resident_number){
    std::vector<int>::iterator ref;
    // Find resident. 
    ref = std::find(current_residents.begin(),current_residents.end(),resident_number);
    if(ref==current_residents.end()){
        std::cout << "Individual is not found!!! " << std::endl;
    }else{
        // Remove the resident.
        current_residents.erase(ref);
    }
}

void assemble_age_matrix(const std::vector<individual> & residents, std::vector<std::vector<int>> & age_matrix){

    // Assign residents to their age matrix.
    for(int i = 0; i < (int) residents.size(); i++){
        int bracket = residents[i].age_bracket;
        age_matrix[bracket].push_back(i);
    }
}

void assemble_age_matrix(const std::vector<individual> & residents, std::vector<std::vector<std::vector<int>>> & age_matrix){
    
    // Assign residents to their age matrix.
    for(int i = 0; i < (int) residents.size(); i++){
        int bracket = residents[i].age_bracket;
        int com_ref = residents[i].current_community;
        age_matrix[com_ref][bracket].push_back(i);
    }
}

// AGE MATRIX INFORMATION
void remove_resident_age_matrix(int resident_number, const individual & resident , std::vector<std::vector<int>> & age_matrix){
    // Find where resident is in the age_matrix.
    int bracket = resident.age_bracket;
    std::vector<int>::iterator ref;
    ref = std::find(age_matrix[bracket].begin(),age_matrix[bracket].end(),resident_number);
        // Was the resident found succesfully... 
        if(ref==age_matrix[bracket].end()){
            std::cout << "Individual is not found!!! " << std::endl;
        }else{
        // Remove the individual from the age_matrix.
            age_matrix[bracket].erase(ref);
        }
        
}


void add_resident_age_matrix(int resident_number, const individual & resident, std::vector<std::vector<int>> & age_matrix){
    // Find where resident is in the age_matrix.
    int bracket = resident.age_bracket;
    age_matrix[bracket].push_back(resident_number);
}


int calculate_SEIR(std::vector<int> &S, std::vector<int> &E,std::vector<int> &I, std::vector<int> &R, const std::vector<individual> & residents){
   int total_infected = 0;
    // This function returns all infected and recovered in each age_bracket.
    for(int i = 0; i < (int) residents.size(); i++){
        const individual & person = residents[i];
            if(person.covid.infection_status=='S')
            {
                int age_bracket = person.age_bracket;
                S[age_bracket] += 1;

            }else if(person.covid.infection_status=='I')
            {
                int age_bracket = person.age_bracket;
                total_infected += 1;
                I[age_bracket] += 1;
            }else if(person.covid.infection_status=='R')
            {
                int age_bracket = person.age_bracket;
                R[age_bracket] +=1;

            }else if(person.covid.infection_status=='E')
            {
                total_infected += 1;
                int age_bracket = person.age_bracket;
                E[age_bracket] += 1;
            }   
    }
    return total_infected;

}

int calculate_SEIR(std::vector<std::vector<int>> &S, std::vector<std::vector<int>> &E,std::vector<std::vector<int>> &I, std::vector<std::vector<int>> &R, const std::vector<individual> & residents){
     int total_infected = 0;
    // This function returns all infected and recovered in each age_bracket.
    for(int i = 0; i < (int) residents.size(); i++){
        const individual & person = residents[i];
        int age_bracket = person.age_bracket;
        // int community_num = person.current_community; // Or community id? 
        int community_num = person.community_id;
            if(person.covid.infection_status=='S')
            {  
                S[community_num][age_bracket] += 1;
            }else if(person.covid.infection_status=='I')
            {   
                total_infected += 1;
                I[community_num][age_bracket] += 1;
            }else if(person.covid.infection_status=='R')
            {   
                R[community_num][age_bracket] +=1;
            }else if(person.covid.infection_status=='E')
            {
                total_infected += 1;
                E[community_num][age_bracket] += 1;
            }   
    }
    return total_infected;
}

int calculate_SEIR_vaccinated(std::vector<int> &S, std::vector<int> &E,std::vector<int> &I, std::vector<int> &R, 
        std::vector<int> &S_V, std::vector<int> &E_V,std::vector<int> &I_V, std::vector<int> &R_V,std::vector<int> &S_V2, std::vector<int> &E_V2,std::vector<int> &I_V2, std::vector<int> &R_V2,const std::vector<individual> & residents){
   int total_infected = 0;
    // This function returns all infected and recovered in each age_bracket.
    for(individual person: residents){
        vaccine_type type = person.vaccine_status.get_type();
        if(type==vaccine_type::none){
            if(person.covid.infection_status=='S')
            {
                int age_bracket = person.age_bracket;
                S[age_bracket] += 1;

            }else if(person.covid.infection_status=='I')
            {
                int age_bracket = person.age_bracket;
                total_infected += 1;
                I[age_bracket] += 1;
            }else if(person.covid.infection_status=='R')
            {
                int age_bracket = person.age_bracket;
                R[age_bracket] +=1;

            }else if(person.covid.infection_status=='E')
            {
                total_infected += 1;
                int age_bracket = person.age_bracket;
                E[age_bracket] += 1;
            }   
        }else if(person.vaccine_status.get_dose()==1){
                if(person.covid.infection_status=='S')
            {
                int age_bracket = person.age_bracket;
                S_V[age_bracket] += 1;

            }else if(person.covid.infection_status=='I')
            {
                int age_bracket = person.age_bracket;
                total_infected += 1;
                I_V[age_bracket] += 1;
            }else if(person.covid.infection_status=='R')
            {
                int age_bracket = person.age_bracket;
                R_V[age_bracket] +=1;

            }else if(person.covid.infection_status=='E')
            {
                total_infected += 1;
                int age_bracket = person.age_bracket;
                E_V[age_bracket] += 1;
            }  
        }else{
            if(person.covid.infection_status=='S')
        {
            int age_bracket = person.age_bracket;
            S_V2[age_bracket] += 1;

        }else if(person.covid.infection_status=='I')
        {
            int age_bracket = person.age_bracket;
            total_infected += 1;
            I_V2[age_bracket] += 1;
        }else if(person.covid.infection_status=='R')
        {
            int age_bracket = person.age_bracket;
            R_V2[age_bracket] +=1;

        }else if(person.covid.infection_status=='E')
        {
            total_infected += 1;
            int age_bracket = person.age_bracket;
            E_V2[age_bracket] += 1;
            
        }
        }
    }
    return total_infected;

}



int calculate_total(std::vector<individual> &residents, std::vector<household> &houses, char ref){

    int val = 0;
    for(int i = 0; i < (int) houses.size(); i++){
        for(int j=0; j < (int) houses[i].current_residents.size();j++){
            int res_number = houses[i].current_residents[j];
            if(residents[res_number].covid.infection_status==ref){
                val+=1;
            }
        }
    }
    return val;
}

int calculate_total(std::vector<individual> &residents, char ref){

    int val = 0;
    for(int i = 0; i < (int) residents.size(); i++){
        // Where are they. 
        individual & person = residents[i];
            if(person.covid.infection_status==ref){
                val+=1;
            }
        
    }
    return val;
}

std::vector<int> calculate_totals(std::vector<individual> &residents, char ref, int num_communities){
    // This function returns all infected in each community.
    std::vector<int> val(num_communities,0);
    for(int i = 0; i < (int) residents.size(); i++){
        individual & person = residents[i];
            if(person.covid.infection_status==ref){
                int c_num = person.community_id;
                val[c_num] += 1;
            }
        
    }
    return val;
}

