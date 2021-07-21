#include "ibm_simulation.h"
#include "individual.h"
#include "households.h"
#include <algorithm>
#include "nbinrnd.h"

static std::random_device rd;
static std::default_random_engine generator(rd());
static std::uniform_real_distribution<double> genunf_std(0.0,1.0);

// // Constructor for disease model, vector beta's
disease_model::disease_model(std::vector<double> beta_H_in, std::vector<double> beta_C_in, std::vector<std::vector<double>> contact_matrix_in, std::vector<double> b,std::vector<double> w):beta(beta_H_in),beta_C(beta_C_in){
    
    contact_matrix = contact_matrix_in; 
    gen_tau_E = std::exponential_distribution<double>(1.0/5.0);
    gen_tau_R = std::exponential_distribution<double>(1.0/2.5);
    gen_tau_S = std::exponential_distribution<double>(1.0/2.5);
    gen_tau_isolation = std::piecewise_constant_distribution<double>(b.begin(),b.end(),w.begin()); // When are they isolated.
    
};

// // Constructor for disease model, scalar beta's
disease_model::disease_model(double beta_H_in,double beta_C_in, std::vector<std::vector<double>> contact_matrix_in,std::vector<double> b,std::vector<double> w):beta(std::vector<double>(contact_matrix_in.size(),beta_H_in)),beta_C(std::vector<double>(contact_matrix_in.size(),beta_C_in)){
    
    contact_matrix = contact_matrix_in; 
    gen_tau_E = std::exponential_distribution<double>(1.0/5.0);
    gen_tau_R = std::exponential_distribution<double>(1.0/2.5);
    gen_tau_S = std::exponential_distribution<double>(1.0/2.5);
    gen_tau_isolation = std::piecewise_constant_distribution<double>(b.begin(),b.end(),w.begin()); // When are they isolated.

};

#ifdef NEW_COVID_MODEL
//  Covid model Age stratified individual contacts. (ASCM - age stratified contact model)
double disease_model::covid_ascm(std::vector<individual> & residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double t0, double t1, double dt, std::vector<size_t> & E, std::vector<size_t> & I,std::vector<size_t> & newly_symptomatic){
//  Evaluates a covid model based on individual contacts from t0 to t1 with timestep of dt.
//  Returns the current time value.

    // Error check the size of the age_ref vs the contact matrix. Throw error if not the same size.
    if(contact_matrix.size()!=age_ref.size()){
        throw std::logic_error("The contact matrix and age reference matrix are not the same size. Check the dimensions of your stratification!");
    }

    if(t0 >= t1){
        throw std::logic_error("Incorrect time specified. t_0 is larger than (or equal to) t_1 in disease_model::covid_ascm.");
    }
    
double t = t0;
// Timestep until the next step would be equal or pass the end time.
    while(t < t1){
        if(t+dt > t1){
            //  Complete the last time step, ensuring that it finishes exactly at tend.
            dt = t1 - t;
        }
        // Run one step of the simulation.
        t = covid_one_step_ascm(residents, houses, age_ref, t, dt, E, I, newly_symptomatic);
        
    }
// Return the time value, should always be tend.
return t;
}

 
bool disease_model::distribution_infected_update(individual & person, size_t & ind_number, std::vector<size_t> & newly_symptomatic, double & t){
    
    // This function is in charge of determining the updates of infected individuals.
    bool symptom_latch = person.covid.check_symptoms; // Define a symptom latch to control that individuals develop symptoms once.
    bool develop_symptoms = (symptom_latch)&&(t > person.covid.time_of_symptom_onset)&&(!person.covid.asymptomatic); // Will they develop symptoms.
    bool recovery = (t > person.covid.time_of_recovery);
    
    if(recovery){
        recover_individual(person); // Do not require t in this function as we can assign time of recovery at point of infection.
    }else if(develop_symptoms){
        // Determine if they become symptomatic (also severity?) // Can be determined at infection time.
        newly_symptomatic.push_back(ind_number); // The individuals that have developed symptoms today.
        person.covid.check_symptoms = false; // Activate latch for this individual.
    } // I am giving precedence to recovering here. If you recover before you develop symptoms, because dt is too large, thats what happens (super unlikely but).

    return recovery;
}


bool disease_model::distribution_exposed_update(individual & person, size_t & ind_number, std::vector<size_t> & newly_infected, double & t){
    // This function returns true if you move from exposed to infected.
    bool infected = (t > person.covid.time_of_infection);
        if(infected){
            infect_individual(person);
            newly_infected.push_back(ind_number); // Must pass in individual number here.
        }
    return infected;
}

//#define NEW_INFECTION_ASCM
#ifdef NEW_INFECTION_ASCM
void disease_model::infection_ascm(double t, individual &infected_individual, household & house, std::vector<individual> & residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double dt, std::vector<size_t> & newly_exposed){
    // We will get the contacts and then infect them.

    // Beta is the probability of infection given contact.
    // Assume they have contact with everyone in the household.
    // Incorporate all of the household information as well for the quarantine.
    
    // This function will alter the infection status of any contacts.
    int     individual_bracket = infected_individual.age_bracket; // Infected individuals age bracket.
    double  ind_transmission = infected_individual.covid.transmissibility*beta[individual_bracket];
    double  ind_community_transmission = (!(infected_individual.isolated||infected_individual.quarantined))*infected_individual.covid.transmissibility*beta[individual_bracket];
    (void) ind_community_transmission;
    
    for(int & house_contact:house.current_residents){
        individual & contact  = residents[house_contact];
        double r = genunf_std(generator);
        double prob_transmission = ind_transmission*contact.vaccine_status.get_susceptibility(t); // Do we multiply by a check on isolation?
        
        if(r < prob_transmission*dt){ // Hack yeah (this part needs the dt to account for the timestep contact, household contacts are all day, whereas other contacts might not be)
                infected_individual.who_infected.push_back(house_contact);
                newly_exposed.push_back(house_contact);
                expose_individual(contact,t);
        };
    }
    
    // Get contacts. (Could skip if isolated etc?)
    std::vector<int> community_contacts; community_contacts.reserve(20); // Magic number to help with reallocation.
    
    for(int & contact_ref: community_contacts){
        individual & contact = residents[contact_ref]; // Define which community member was contacted.
        bool house_q = houses[contact.current_home].external_quarantine || houses[contact.current_home].quarantined;
        
        // If you are isolated, or are in a house that is quarantined, you cant be accessed. This has to be incorporated.
        // I need to check who is in the community. If you are quarantined you cant be accessed. Should this be done by changing age_ref?  Or incorpoated here. Probably age_ref.
        // Can only infect susceptible contacts
        if(contact.isolated || contact.external_isolation|| house_q){

        }else{
        // Track all contacts (We can then implement how good contact tracing is currently disabled)

        if(contact.covid.infection_status == 'S'){

        // Do they get infected.
        double r = genunf_std(generator);
        double prob_transmission = infected_individual.covid.transmissibility*beta_C[individual_bracket]*contact.vaccine_status.get_susceptibility(t);
            if(r < prob_transmission){   // This part doesnt need the dt, because the dt is taken into account earlier by limiting the number of contacts per timesteps.
                newly_exposed.push_back(contact_ref);
                expose_individual(contact,t);
                infected_individual.who_infected.push_back(contact_ref);
                contact.covid.cluster_number = cluster_number;
            }
            
        }

    }
    }
}

#else

void disease_model::infection_ascm(double t, individual &infected_individual, household & house, std::vector<individual> & residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double dt, std::vector<size_t> & newly_exposed){
    
    // It is definitely worth splitting this function up into generating contacts and infections.
    
    // This function does not change anything about house or the infected individual.
    // Beta is the probability of infection given contact.
    // Assume they have contact with everyone in the household.
    // Incorporate all of the household information as well for the quarantine.
    
    // This function will alter the infection status of any contacts. 
    int individual_bracket = infected_individual.age_bracket; // Infected individuals age bracket.
    int cluster_number = infected_individual.covid.cluster_number;
    // External isolation or quarantine.
    if(infected_individual.external_isolation||infected_individual.external_quarantine||house.external_quarantine){
        // No interactions
    }else{
    // Household transmission.
    // The dt has been included here as it means that the probability doesnt change when time step changes. Its a hack that works perfectly.
    // We do not check if isolated here as they are in the house. They will be interacting.
    for(size_t i = 0; i < house.current_residents.size(); i++){

        individual & contact  = residents[house.current_residents[i]];
        
        int     contact_ref = house.current_residents[i];

        // If contact is quarantined or isolated they shouldnt be contacted. This is the household, so should be external isolation and quarantine.
        if(contact.external_isolation||contact.external_quarantine){
            // Note that we do not have to check if the contacts house is quarantined because this is a household contact (they live in the same place)
        }   else if(contact.covid.infection_status == 'S'){
            // Can only infect susceptible contacts
            // Do they get infected.
            double r = genunf_std(generator);
            double prob_transmission = infected_individual.covid.transmissibility*beta[individual_bracket]*contact.vaccine_status.get_susceptibility(t);
            if(std::isnan(prob_transmission)){
                throw std::logic_error("Probability of infection is Nan.");
            }
            if(r < prob_transmission*dt){ // Hack yeah (this part needs the dt to account for the timestep contact, household contacts are all day, whereas other contacts might not be)
                    infected_individual.who_infected.push_back(contact_ref);
                    newly_exposed.push_back(contact_ref);
                    expose_individual(contact,t);
                contact.covid.cluster_number = cluster_number;
            };
        }
    };

    //  Community transmission. The person will contact people depending upon the contact matrix. For each component of the contact matrix we sample how many contacts they make from each age bracket.
    //  This checks their current house. But what if theyre not in their current house. Could probably still interact. A
    if(infected_individual.isolated || infected_individual.external_isolation || infected_individual.external_quarantine|| house.quarantined){
        // Do nothing.

    }else{
        
    for(size_t age_strata = 0; age_strata < contact_matrix.size(); age_strata++){
        
        // Average number of contacts per day. 
        double mu_contacts = contact_matrix[individual_bracket][age_strata];

        // Note that contact matrix is square so doesnt matter which dimension has their size checked.
        std::poisson_distribution<int> gen_num_contacts(mu_contacts*dt); // Create the probability distributon for the number of contacts. 

        // The dt scale is so that it makes sense.
        size_t num_in_strata = age_ref[age_strata].size(); // Must be obtained from the reference to the age matrix.

        if(num_in_strata == 0){
            continue;    // Skip strata as no-one is in it.
        }

        // I dont bother checking if the individual interacts with themselves (This is very unlikely to happen)
        int number_comm_contacts = gen_num_contacts(generator);
            for(int i = 0; i < number_comm_contacts; i++){
                // Should replace rand with appropriate distribution from <random>
                int contact_ref = age_ref[age_strata][rand()%num_in_strata]; // Sample with replacement who is contacted from the community, could be switched to without replacement(surely wont matter).
                
                // This is where we can finish generating contacts, the rest of the loop calculated whether they are infected. Should this be split ?
                individual & contact = residents[contact_ref]; // Define which community member was contacted.
                bool house_q = houses[contact.current_home].external_quarantine || houses[contact.current_home].quarantined;
                // If you are isolated, or are in a house that is quarantined, you cant be accessed. This has to be incorporated.
                // I need to check who is in the community. If you are quarantined you cant be accessed. Should this be done by changing age_ref?  Or incorpoated here. Probably age_ref.
                // Can only infect susceptible contacts
                if(contact.isolated || contact.external_isolation|| house_q){

                }else{
                // Track all contacts (We can then implement how good contact tracing is currently disabled)

                if(contact.covid.infection_status == 'S'){

                // Do they get infected.
                double r = genunf_std(generator);
                double prob_transmission = infected_individual.covid.transmissibility*beta_C[individual_bracket]*contact.vaccine_status.get_susceptibility(t);
                    if(r < prob_transmission){   // This part doesnt need the dt, because the dt is taken into account earlier by limiting the number of contacts per timesteps.
                        newly_exposed.push_back(contact_ref);
                        expose_individual(contact,t);
                        infected_individual.who_infected.push_back(contact_ref);
                        contact.covid.cluster_number = cluster_number;
                    }
                    
                }
                }
            }
    }
    }
    }

}
#endif

double  disease_model::covid_one_step_ascm(std::vector<individual>& residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double t0, double dt, std::vector<size_t> & E, std::vector<size_t> & I,std::vector<size_t> & newly_symptomatic){
    // The ordering of this function is important. Infected individuals, then exposed.
    // There is a shift in memory location for each true statement that comes from std::remove_if. This is faster than creating a new vector (checked with example).
    
    // Reserve memory for new infections... Small but should save from reallocation costs which can be significant.
    std::vector<size_t> newly_exposed; newly_exposed.reserve(3000); // 3000 is a magic number... could try different values. Somehow relate to dt could be possible.
    std::vector<size_t> newly_infected; newly_infected.reserve(3000); // Could make dependent upon the number of exposed individuals.
    
    // Loop over all infected individuals, removing those that recover. This will alter a vector newly_exposed to have the newly exposed individuals.
    auto recovered_it = std::remove_if(I.begin(),I.end(),[&](size_t & ind_ref)->bool{
       
        // Infected individual
        individual & person = residents[ind_ref];
        
            //Error check.
            if(person.covid.infection_status!='I'){
                throw std::logic_error("Individual in I vector does not match infection status.");
            }
            
        // Infected individuals household.
        household  & house  = houses[person.current_home];
        
        // Household infection model (currently a part of infection_ascm).
        
        // Community infection model.
        infection_ascm(t0, person, house, residents, houses, age_ref, dt, newly_exposed);
        
        return distribution_infected_update(person, ind_ref, newly_symptomatic, t0); // If infected, do they recover, this function must alter the individual disease status.
    });
    
    // Loop over Exposed individuals. Note that this does not include newly_exposed individuals.
    auto infected_it = std::remove_if(E.begin(),E.end(),[&](size_t & ind_ref)->bool{
        
        // Get person of interest.
        individual & person = residents[ind_ref];
        
            //Error check.
            if(person.covid.infection_status!='E'){
                throw std::logic_error("Individual in I vector does not match infection status.");
            }
        
        bool infected = distribution_exposed_update(person,ind_ref,newly_infected, t0); // Distribution exposed update should now update individuals to be infected., this function must alter the individual disease status.
        
        return infected;
        
    });
    
    // We have a vector of newly exposed individuals and an iterator, recovered_it, that points to the location of recovered individuals.
    I.erase(recovered_it,I.end()); // Remove the recovered individuals from the vector of infections (could be done at the end).
    I.insert(I.end(),newly_infected.begin(),newly_infected.end()); // Adds the removed end of E to the end of I.
    
    // We cannot insert the remove_if values to the end of I as the value after remove_if is unspecified.
    E.erase(infected_it,E.end()); // Remove the exposed individuals that are now infected.
    E.insert(E.end(),newly_exposed.begin(),newly_exposed.end()); // Insert the newly exposed individuals to the end of exposed individuals.
    
    // The probability of tranisitioning into a new compartment is now dependent upon the order of the individuals in E and I. This should be better for branch prediction. It is still stochastic so its not just perfect order, but its better! (Very likely to be false at the end).

    return t0 + dt;
}
#endif

[[deprecated("disease_model::seed_exposure(individual ) is deprecated, please use the disease_model::seed_exposure(individual , double )")]]
void disease_model::seed_exposure(individual &resident){
    expose_individual(resident);
}

void disease_model::seed_exposure(individual & resident, double & t){
    expose_individual(resident, t);
}

[[deprecated("disease_model::expose_individual(individual ) is deprecated, please use the disease_model::expose_individual(individual , double )")]]
inline void disease_model::expose_individual(individual & resident){
    resident.covid.infection_status = 'E';
}
//[[deprecated("disease_model::infect_individual(individual ) is deprecated, please use the disease_model::infect_individual(individual , double )")]]
inline void disease_model::infect_individual(individual & resident){
    resident.covid.infection_status = 'I';
}
//[[deprecated("disease_model::recover_individual(individual ) is deprecated, please use the disease_model::recover_individual(individual , double )")]]
inline void disease_model::recover_individual(individual & resident){
    resident.covid.infection_status = 'R';
}
//[[deprecated("disease_model::susceptible_individual(individual ) is deprecated, please use the disease_model::susceptible_individual(individual , double )")]]
inline void disease_model::susceptible_individual(individual & resident){
    resident.covid.infection_status = 'S';
}

inline void disease_model::expose_individual(individual & resident, double & t){
    
    resident.covid.infection_status = 'E';
    resident.covid.time_of_exposure = t;
    resident.covid.time_of_infection = resident.covid.time_of_exposure + gen_tau_E(generator);
    resident.covid.time_of_symptom_onset = resident.covid.time_of_infection + gen_tau_S(generator);
    resident.covid.time_of_recovery = resident.covid.time_of_symptom_onset + gen_tau_R(generator);
    
    resident.time_isolated = resident.covid.time_of_symptom_onset - gen_tau_isolation(generator);
    
    // Determine if the individual will be asymptomatic and the severity of the disease.
    std::uniform_real_distribution<double> gen_r(0.0,1.0);
    double r = gen_r(generator);
    double prob_asymptomatic = resident.vaccine_status.get_probability_asymptomatic(t);
    bool   asymptomatic = r > prob_asymptomatic; // To fix the asymptomatic should just be symptomatic. It was defined the opposite way. woops.
    resident.covid.asymptomatic = asymptomatic;
    
    // Severity check (Currently in model of care).
    resident.covid.severe = false; // This is currently disabled.
    
    // Determine the transmissibility of the individual.
    resident.covid.transmissibility = resident.vaccine_status.get_transmissibility(t,asymptomatic?1:0); // This will remain over the course of their infection.
    
    // Set statistics for tracking.
    resident.infection_statistics.set_infection_stats(resident.vaccine_status.get_type(), resident.vaccine_status.get_dose(), resident.vaccine_status.get_time_of_vaccination());
    
}


// This code is used to check the R0.
double disease_model::covid_ascm_R0(std::vector<individual> & residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double t0, double t1, double dt, std::vector<size_t> & E, std::vector<size_t> & I,std::vector<size_t> & newly_symptomatic){
//  Evaluates a covid model based on individual contacts from t0 to t1 with timestep of dt.
//  Returns the current time value.

    // Error check the size of the age_ref vs the contact matrix. Throw error if not the same size.
    if(contact_matrix.size()!=age_ref.size()){
        throw std::logic_error("The contact matrix and age reference matrix are not the same size. Check the dimensions of your stratification!");
    }

    if(t0 >= t1){
        throw std::logic_error("Incorrect time specified. t_0 is larger than (or equal to) t_1 in disease_model::covid_ascm.");
    }
    
double t = t0;
// Timestep until the next step would be equal or pass the end time.
    while(t < t1){
        if(t+dt > t1){
            //  Complete the last time step, ensuring that it finishes exactly at tend.
            dt = t1 - t;
        }
        // Run one step of the simulation.
        t = covid_one_step_ascm_R0(residents, houses, age_ref, t, dt, E, I, newly_symptomatic);
        
    }
// Return the time value, should always be tend.
return t;
}

double  disease_model::covid_one_step_ascm_R0(std::vector<individual>& residents, std::vector<household> & houses, std::vector<std::vector<int>> & age_ref, double t0, double dt, std::vector<size_t> & E, std::vector<size_t> & I,std::vector<size_t> & newly_symptomatic){
    // E and I are all that is needed to generate incidence? (Do I want to add E and I to the disease_model class?)

    // The ordering of this function is important. Infected individuals, then exposed.
    // There is a shift in memory location for each true statement that comes from std::remove_if. This is faster than creating a new vector (checked with example).
    // Reserve memory for new infections... Small but should save from reallocation costs which can be significant.
    std::vector<size_t> newly_exposed; newly_exposed.reserve(3000); // 3000 is a magic number... could try different values. Somehow relate to dt could be possible.
    std::vector<size_t> newly_infected; newly_infected.reserve(1); // Could make dependent upon the number of exposed individuals.
    
    // Loop over all infected individuals, removing those that recover. This will alter a vector newly_exposed to have the newly exposed individuals.
    auto recovered_it = std::remove_if(I.begin(),I.end(),[&](size_t & ind_ref)->bool{
       
        // Infected individual
        individual & person = residents[ind_ref];
        
            //Error check.
            if(person.covid.infection_status!='I'){
                throw std::logic_error("Individual in I vector does not match infection status.");
            }
            
        // Infected individuals household.
        household  & house  = houses[person.current_home];
        
        // Household infection model (currently a part of infection_ascm).
        
        // Community infection model.
        infection_ascm(t0, person, house, residents, houses, age_ref, dt, newly_exposed);
        
        return distribution_infected_update(person, ind_ref, newly_symptomatic, t0); // If infected, do they recover, this function must alter the individual disease status.
    });
    
    // Loop over Exposed individuals. Note that this does not include newly_exposed individuals.
    auto infected_it = std::remove_if(E.begin(),E.end(),[&](size_t & ind_ref)->bool{
        
        // Get person of interest.
        individual & person = residents[ind_ref];
        
            //Error check.
            if(person.covid.infection_status!='E'){
                throw std::logic_error("Individual in I vector does not match infection status.");
            }
        
        bool infected = distribution_exposed_update(person,ind_ref,newly_infected, t0); // Distribution exposed update should now update individuals to be infected., this function must alter the individual disease status.
        
        return infected;
        
    });
    
    // We have a vector of newly exposed individuals and an iterator, recovered_it, that points to the location of recovered individuals.
    I.erase(recovered_it,I.end()); // Remove the recovered individuals from the vector of infections (could be done at the end).
    I.insert(I.end(),newly_infected.begin(),newly_infected.end()); // Adds the removed end of E to the end of I.
    
    // We cannot insert the remove_if values to the end of I as the value after remove_if is unspecified.
    E.erase(infected_it,E.end()); // Remove the exposed individuals that are now infected.
    // Do not add newly exposed individuals, we just want them to be exposed.
    
    return t0 + dt;
}
