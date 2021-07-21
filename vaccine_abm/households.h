# ifndef HOUSEHOLDS_H
# define HOUSEHOLDS_H
# include <vector>
# include <iostream>
# include <random>

class individual;   // Class definition required for individuals.

class household
{
public:
    //  Constructor functions for households. 
    household();
    household(int , int , int , std::vector<individual> &,const std::vector<double> &,const std::vector<double> pi,vaccine_parameters &);
    
    //      Member Data
    int     household_size;         //  This is currently not required, but why not.
    bool    quarantined = false;    //  Is the household quarantined?
    bool    external_quarantine = false; // Are they externally quarantined.
    double  time_quarantined;       //  How long has the house been quarantined.    
    std::vector<int> current_residents; // Who is currently in the house ? 

    // Member functions.
    void add_resident(int);         //  Add resident to the house.
    void remove_resident(int);

};

class community
{
public:
    //  Constructor functions.
    community();
    community(int comm_id, int n_houses ,double house_average, std::vector<individual> & residents ,std::vector<household> & houses, const std::vector<double> & age_brackets,const std::vector<double> pi, vaccine_parameters &);     //  Constructor function for the community
    //  Member data.
    int     population_size = 0;            // Total number of individuals.
    int     cumulative_incidence = 0;       // Cumulative incidence. 

    bool    quarantined = false;            // Quarantine the community. 
    double  time_quarantined;               //  How long has the house been quarantined.    

};

int calculate_total_susceptible(std::vector<household> &);
int calculate_total_exposed(std::vector<household> &);
int calculate_total_infection(std::vector<household> &);
int calculate_total_recovered(std::vector<household> &);
void assemble_age_matrix(const std::vector<individual> & residents, std::vector<std::vector<std::vector<int>>> & age_matrix);
void assemble_age_matrix(const std::vector<individual> & residents, std::vector<std::vector<int>> & age_matrix);

int calculate_total(std::vector<individual> &residents, std::vector<household> &houses,char);
int calculate_total(std::vector<individual> &residents,char);
std::vector<int> calculate_totals(std::vector<individual> &residents, char ref, int num_communities);
int calculate_SEIR(std::vector<int> &S, std::vector<int> &E,std::vector<int> &I, std::vector<int> &R, const std::vector<individual> & residents);
int calculate_SEIR(std::vector<std::vector<int>> &S, std::vector<std::vector<int>> &E,std::vector<std::vector<int>> &I, std::vector<std::vector<int>> &R, const std::vector<individual> & residents);

// Age matrix information
void remove_resident_age_matrix(int resident_number,const individual &, std::vector<std::vector<int>> & age_matrix);
void add_resident_age_matrix(int resident_number, const individual & resident, std::vector<std::vector<int>> & age_matrix);
void update_age_matrix(int, int, std::vector<std::vector<int>> &, std::vector<std::vector<int>> &);
int calculate_SEIR_vaccinated(std::vector<int> &S, std::vector<int> &E,std::vector<int> &I, std::vector<int> &R, std::vector<int> &S_V, std::vector<int> &E_V,std::vector<int> &I_V, std::vector<int> &R_V,std::vector<int> &, std::vector<int> &,std::vector<int> &, std::vector<int> &,const std::vector<individual> & residents);
# endif
