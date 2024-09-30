////
////  main.cpp
////  Experiment Genetic Algorithm
////
////  Created by Angel Cruz Mora on 27/08/24.
////

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <numeric>
#include <algorithm>
#include <tuple>
#include <time.h>
#include <random>

#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>


// Configuration file path
std::string CONFIGURATION_PATH = "/Users/angelcruz/Documents/C++/Experimento - AG/Experimento - AG/config.txt";


// Parameters Configuration
std::unordered_map<std::string, float> PARAMETERS =
{
    {"POPULATION", 100},
    {"GENES", 20},
    {"MAX_ITERATIONS", 50},
    {"RUNS_NUMBER", 30},
    {"SELECTION", 0},
    {"TOURNAMENT_INDIVIDUALS", 10},
    {"TRUNCATION_RATE", 0.5},
    {"SELECTION_RATE", 0.5},
    {"CROSSOVER", 0},
    {"CROSSOVER_PROBABILITY", 0.5},
    {"MUTATION", 0},
    {"MUTATION_PROBABILITY", 0.9},
    {"ERROR_FLAG", 0}
};

// Auxiliary methdos to generate random values.
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);


/*  AUXILIARY FUNCTIONS  */
float computeMean(std::vector<int> values);
void load_configuration();
auto init_parameters() -> std::vector<std::vector<int>>;
auto compute_fitness(std::vector<std::vector<int>> population) -> std::vector<int>;
auto sort_individuals(std::vector<std::vector<int>> population, std::vector<int> &fitness) -> std::vector<std::vector<int>>;
/* SELECTION METHODS */
auto roulette_selection(std::vector<std::vector<int>> population, std::vector<int> fitness) -> std::vector<std::vector<int>>;
auto tournamentSelection(std::vector<std::vector<int>> population, std::vector<int> fitness) -> std::vector<std::vector<int>>;
auto universalStochasticSelection(std::vector<std::vector<int>> population, std::vector<int> fitness) -> std::vector<std::vector<int>>;
auto selectionByTruncation(std::vector<std::vector<int>> population, std::vector<int> fitness) -> std::vector<std::vector<int>>;
/* CROSSOVER METHODS */
auto two_points_crossover(std::vector<int> ind1, std::vector<int> ind2) -> std::vector<std::vector<int>>;
auto uniform_crosssover(std::vector<int> ind1, std::vector<int> ind2) -> std::vector<std::vector<int>>;
auto randomMask(std::vector<int> ind1, std::vector<int> ind2) -> std::vector<std::vector<int>>;
/* MUTATION METHODS */
void FlipBits(std::vector<std::vector<int>> individuals);
void tripleMutation(std::vector<std::vector<int>> individuals);
void reArrangeBits(std::vector<std::vector<int>> individuals);


/*------------------------------------------------------------------------------------------------------------*/


auto oneMaxExperiment(std::vector<int> best_ind, float &best_fitness)
{
    int ite = 0;
    
    // Define the number of runs
    int MAX_ITERATIONS = PARAMETERS["MAX_ITERATIONS"];
    // Number of Genes
    int GENES = PARAMETERS["GENES"];
    
    // Defininig the operators of the algorithm
    int selectionMethod = PARAMETERS["SELECTION"];
    int crossoverMethod = PARAMETERS["CROSSOVER"];
    int mutationMethod = PARAMETERS["MUTATION"];
    
    
    // Define a vector of vectors to store the initial population
    std::vector<std::vector<int>> population = init_parameters();
    // Define a vector to store the fitness of each individual
    std::vector<int> fitness;
    // Define a vector of vectors to store the next generation
    std::vector<std::vector<int>> offspring;
    // Define a vector of vectors to store the parents of the next individuals
    std::vector<std::vector<int>> parents;
    // Define a vector of vectors to store the next individuals
    std::vector<std::vector<int>> aux;
    
    
    // Define an int to use it as flag, since we generate two individuals at a time we'll have to drop one
    int odd_pop = (int)PARAMETERS["POPULATION"] % 2;
    
    // Main Loop
    while(true)
    {
        // Compute fitness of each individual
        fitness = compute_fitness(population);
        
        // Sort individuals based on their fitness
        population = sort_individuals(population, fitness);
        
        // If statement to get off the while loop
        if (ite == MAX_ITERATIONS || fitness[0] == GENES)
        {
            break;
        }
        
        // For loop to generate the next generation
        for(int i = 0; i < PARAMETERS["POPULATION"] / 2 + odd_pop; i++)
        {
            // Selection Procedure
            switch (selectionMethod) {
                case 0:
                    parents = roulette_selection(population, fitness);
                    break;
                case 1:
                    parents = tournamentSelection(population, fitness);
                    break;
                case 2:
                    parents = universalStochasticSelection(population, fitness);
                    break;
                case 3:
                    parents = selectionByTruncation(population, fitness);
                    
                default:
                    parents = roulette_selection(population, fitness);
            }
            
            // Crossover Procedure
            switch (crossoverMethod) {
                case 0:
                    aux = two_points_crossover(parents[0], parents[1]);
                    break;
                case 1:
                    aux = uniform_crosssover(parents[0], parents[1]);
                    break;
                case 2:
                    aux = randomMask(parents[0], parents[1]);
                    break;
                    
                default:
                    aux = two_points_crossover(parents[0], parents[1]);
            }
            
            
            // Mutation Procedure
            switch (mutationMethod) {
                case 0:
                    FlipBits(aux);
                    break;
                case 1:
                    tripleMutation(aux);
                    break;
                case 2:
                    reArrangeBits(aux);
                    break;
                    
                default:
                    break;
            }
            
            // Store the new individuals generated
            offspring.push_back(aux[0]);
            offspring.push_back(aux[1]);
        }
        // If population is an odd number remove the
        // extra individual, since we add two at a time
        if(odd_pop)
        {
            offspring.pop_back();
        }
        
        // Replace the old population with the new population
        population = offspring;
        offspring.clear();
        
        ite++;
        std::cout << "-";
    }
    std::cout << "\n";
    
    //best_ind = population[0];
    best_fitness = fitness[0];
    
    return population[0];
}







int main() {
    srand( (unsigned)time( NULL ) );
    
    // Initiliaze Configuration
    load_configuration();
    
    // Verify we load the configuration successfully
    if(PARAMETERS["ERROR_FLAG"])
    {
        std::cout << "\n\tWe couldn't Configure the Algorithm, please add the config.txt file!\n\n";
        return 0;
    }
    
    // Get the number of Generations
    int RUNS_NUMBER = PARAMETERS["RUNS_NUMBER"];
    
    // Define a vector to store the best individual's fitness of the i-th run
    std::vector<int> bestIndividualsFitness;
    
    // Define auxilary vectors to store the best individual and its fitness
    std::vector<int> aux_vect;
    float aux_int;
    
    // Start the algorithm
    for(int i = 0; i < RUNS_NUMBER; i++)
    {
        std::cout << "\n";
        
        aux_vect = oneMaxExperiment(aux_vect, aux_int);
        
        // Print statistics
        std::cout << "RUN: " << i + 1 << "\n Mejor Individuo. \n\tFitness: " << aux_int << " - ";
        
        // Print the best individual
        for(int j = 0; j < PARAMETERS["GENES"]; j++)
        {
            std::cout << aux_vect[j] << " ";
        }
        
        // Store the best fitness
        bestIndividualsFitness.push_back(aux_int);
    }
    
    float genralMean = computeMean(bestIndividualsFitness);
    
    // Print general statistics
    std::cout << "\n\n\n\t * NUMBER OF RUNS. " << RUNS_NUMBER <<  "\n\t ** EXPERIMENT GENERAL MEAN. " << genralMean;
    std::cout << "\n\n" << std::endl;

    return 0;
}



/* ************************************************************************ */
/**

 
            -- AUXILIARY FUNCTIONS --
 
 */

float computeMean(std::vector<int> values){
    
    float mean = (float)accumulate(values.begin(), values.end(), 0) / values.size();
    return mean;
}

void load_configuration()
{
    /*
        Function to load the configuration file and use it to performe the run
    */
    
    // Define variables to extract content from text file
    std::string line;
    std::string token;
    std::string key;
    
    // Open input file
    std::ifstream inputFile(CONFIGURATION_PATH);
    
    // Verify if the file was open
    if (!inputFile.is_open())
    {
        std::cout << "\n\tWe couldn't find the config.txt file!\n";
        PARAMETERS["ERROR_FLAG"] = 1;
        return;
    }

    // Extract the content
    while (std::getline(inputFile, line))
    {
        if (line == "" || line.find("#") == 0)
        {
            continue;
        }

        // Split line by spaces
        std::istringstream iss(line);
        
        // Get Parameter Key
        iss >> token;
        key = token;

        // Iterate along the elements until we get the last element, which is the value
        while (iss >> token) {
            continue;
        }
        
        // Assign Parameter value
        try
        {
            PARAMETERS[key] = std::stof(token);
        }
        catch(const std::invalid_argument& e)
        {
            std::cerr << "Invalid value for " << key << " : " << token << "\n";
        }
    }
    
    inputFile.close();
}

/* ************************************************************************ */

auto init_parameters() -> std::vector<std::vector<int>>
{
    /*
        Function to initialize the population, return individuals filled with random values, either 0 or 1.
    */
    
    // Load the Configuration from the config.txt file
    load_configuration();
    
    // Define population vector
    std::vector<std::vector<int>> pobl;
    std::vector<int> ind;
    
    for(int i = 0; i < PARAMETERS["POPULATION"]; i++){
        for(int j = 0; j < PARAMETERS["GENES"]; j++){
            ind.push_back(rand() % 2);
        }
        pobl.push_back(ind);
        ind.clear();
    }
    return pobl;
}

/* ************************************************************************ */

auto compute_fitness(std::vector<std::vector<int>> population) -> std::vector<int>
{
    /*
        Function to compute the fitness of each individual.
     
            Return:
                    A vector<int> with the fitness of each individual.
     */
    std::vector<int> fitness;
    
    for(int i = 0; i < PARAMETERS["POPULATION"]; i++){
        fitness.push_back(accumulate(population[i].begin(), population[i].end(), 0));
    }
    return fitness;
}


/* ************************************************************************ */


auto sort_individuals(std::vector<std::vector<int>> population, std::vector<int> &fitness) -> std::vector<std::vector<int>> //std::vector<std::vector<int>>
{
    /*
        Algorithm to sort the individuals based on its fitness, from max to min.
                Parameters:
                        vector<vector<int>> population, vector<int> fitness
    */
    
    // Auxiliary Variables
    std::vector<int> ind;
    int aux, j;
    
    for(int t = 0; t < population.size(); t++){
        for(int i = 0; i < population.size() - 1; i++){
            j= i+1;
            ind= population[i];
            aux= fitness[i];
            while(j < population.size() && aux < fitness[j]){
                j++;
            }
            fitness[i] = fitness[j - 1];
            fitness[j - 1] = aux;
            population[i] = population[j - 1];
            population[j - 1] = ind;
        }
    }
    
    return population;
}


/* ************************************************************************ */
/**
 
            -- SELECTION METHODS --
        
 */

auto roulette_selection(std::vector<std::vector<int>> population, std::vector<int> fitness) -> std::vector<std::vector<int>>
{
    /*
            Function to perform the parents selection using the Roulette Method.
                Parameters.
                            vector<vector<int>> population;
                            vector<int> fitness;
     
                Returns.
                            vector<vector<int>> parents
                                                        The parents selected.
    */
    std::vector<std::vector<int>> parents;
    int sum_, u, total_fit= accumulate(fitness.begin(), fitness.end(), 0);
    
    
    // Roulette Method
    for(int i = 0; i < 2; i++){
        
        u= dis(gen) * total_fit;
        sum_ = 0;
        
        for(int j = 0; j < PARAMETERS["POPULATION"]; j++){
            sum_ = sum_ + fitness[j];
            
            if(sum_ > u){
                parents.push_back(population[j]);
                break;
            }
        }
    }
    
    
    return parents;
}


auto selectionByTruncation(std::vector<std::vector<int>> population, std::vector<int> fitness) -> std::vector<std::vector<int>>
{
    /*
            Function to perform the parents selection by truncation.
            We will sort the individuals by fitness and then take P porcentage
            of the total population.
                Parameters.
                            vector<vector<int>> population;
                            vector<int> fitness;
     
                Returns.
                            vector<vector<int>> parents
                                                        The parents selected.
    */
    std::vector<std::vector<int>> parents;
    unsigned int aux = PARAMETERS["TRUNCATION_RATE"] * PARAMETERS["POPULATION"];
    std::uniform_int_distribution<> rand_num(0, aux);
    
    for(int i = 0; i < 2; i++)
    {
        parents.push_back(population[rand_num(gen)]);
    }
    
    return parents;
}

auto universalStochasticSelection(std::vector<std::vector<int>> population, std::vector<int> fitness) -> std::vector<std::vector<int>>
{
    /*
            Function to perform the parents selection using the
            Universal Stochastic Selection. By default we will select just two individuals.
                Parameters.
                            vector<vector<int>> population;
                            vector<int> fitness;
     
                Returns.
                            vector<vector<int>> parents
                                                        The parents selected.
    */
    std::vector<std::vector<int>> parents;
    int sum_, total_fit= accumulate(fitness.begin(), fitness.end(), 0);
    int i, aux;
    
    // Total Fitness  divided by the k individuals to select
    int distance = total_fit / 2;
    
    // Generate the starting point between 0 and the distance
    std::uniform_int_distribution<> rand_num(0, distance);
    int start = rand_num(gen);
    
    // Generate the random points
    int points[] = {0, 0};
    
    for(i = 0; i < 2; i++)
    {
        points[i] = start + i * distance;
    }
    
    // Selection Procedure
    for(i = 0; i < 2; i++)
    {
        aux = 0;
        sum_ = fitness[0];
        while(sum_ < points[i])
        {
            aux++;
            sum_ += fitness[aux];
        }
        parents.push_back(population[aux]);
    }
    
    return parents;
}

auto tournamentSelection(std::vector<std::vector<int>> population, std::vector<int> fitness) -> std::vector<std::vector<int>>
{
    /*
            Function to perform the parents selection using the
            Tournament Selecion.
                Parameters.
                            vector<vector<int>> population;
                            vector<int> fitness;
     
                Returns.
                            vector<vector<int>> parents
                                                        The parents selected.
    */
    std::vector<std::vector<int>> parents;
    std::uniform_int_distribution<> rand_num(0, PARAMETERS["POPULATION"] - 1);
    unsigned int rand_ind, max_ind, max_fitness = 0, k;
    
    for(int i = 0; i < 2; i ++)
    {
        k = PARAMETERS["TOURNAMENT_INDIVIDUALS"];
        
        while(k--)
        {
            rand_ind = rand_num(gen);
            
            if(max_fitness < fitness[rand_ind])
            {
                max_fitness = fitness[rand_ind];
                max_ind = rand_ind;
            }
        }
        parents.push_back(population[max_ind]);
    }
    
    return parents;
}


/* ************************************************************************ */
/**
 
            -- CROSSOVER METHODS --
        
 */

auto two_points_crossover(std::vector<int> ind1, std::vector<int> ind2) -> std::vector<std::vector<int>>
{
    /*
        Function to perform the crossover betweeen two individuals.
                
                Parameters.
                            vector<int> ind1;
                            vector<int> ind2;
     
                Return.
                            vector<vector<int>> offspring
    */
    std::uniform_real_distribution<> gen_cxpoint1(1, PARAMETERS["GENES"]);
    std::uniform_real_distribution<> gen_cxpoint2(1, PARAMETERS["GENES"] - 1);
    std::vector<std::vector<int>> offspring;
    
    int point1, point2, temp;
    
    // Generate random points
    point1 = gen_cxpoint1(gen);
    point2 = gen_cxpoint2(gen);
    
    if(point2 >= point1)
    {
        point2++;
    }
    else
    {
        temp = point1;
        point1 = point2;
        point2 = temp;
    }
    
    // Performe the crossover
    for(int i = point1; i < point2; i++){
        temp = ind1[i];
        ind1[i] = ind2[i];
        ind2[i] = temp;
    }
    
    offspring.push_back(ind1);
    offspring.push_back(ind2);
    
    return offspring;
}


auto uniform_crosssover(std::vector<int> ind1, std::vector<int> ind2) -> std::vector<std::vector<int>>
{
    /*
        Function to perform the crossover betweeen two individuals.
                
                Parameters.
                            vector<int> ind1;
                            vector<int> ind2;
     
                Return.
                            vector<vector<int>> offspring
    */
    
    std::vector<std::vector<int>> offspring;
    int aux;
    
    for(int i = 0; i < PARAMETERS["GENES"]; i++){
        if (dis(gen) < PARAMETERS["CROSSOVER_PROBABILITY"]){
            aux = ind1[i];
            ind1[i] = ind2[i];
            ind2[i] = aux;
        }
    }
    
    offspring.push_back(ind1);
    offspring.push_back(ind2);
    
    return offspring;
}


auto randomMask(std::vector<int> ind1, std::vector<int> ind2) -> std::vector<std::vector<int>>
{
    /*
        Function to perform the crossover Random Mask betweeen two individuals.
                
                Parameters.
                            vector<int> ind1;
                            vector<int> ind2;
     
                Return.
                            vector<vector<int>> offspring
    */
    
    std::vector<std::vector<int>> offspring;
    int i, aux;
    
    // Generate Random mask
    int mask[(int)PARAMETERS["GENES"]];
    std::uniform_int_distribution<> rand_num(0, 1);
    
    for(i = 0; i < PARAMETERS["GENES"]; i++)
    {
        mask[i] = rand_num(gen);
    }
    
    // Selection Procedure
    for(i = 0; i < PARAMETERS["GENES"]; i++)
    {
        aux = ind1[i];
        ind1[i] = (mask[i])? ind2[i] : ind1[i];
        ind2[i] = (mask[i])? ind2[i] : aux;
    }
    
    offspring.push_back(ind1);
    offspring.push_back(ind2);
    
    return offspring;
}


/* ************************************************************************ */
/**
 
            -- MUTATION METHODS --
        
 */

void FlipBits(std::vector<std::vector<int>> individuals)
{
    /*
        Function to performe mutation. It flips the bits of the entire individual.
                
                Parameters.
                                vector<vector<int>> individuals
                                float mutpb  -> mutation probability
    */
    
    for(int idx = 0; idx < individuals.size(); idx++)
    {
        
        for(int j = 0; j < PARAMETERS["GENES"]; j++)
        {
            if(dis(gen) < PARAMETERS["MUTATION_PROBABILITY"])
            {
                individuals[idx][j] = 1 - individuals[idx][j];
            }
        }
    }
}


void tripleMutation(std::vector<std::vector<int>> individuals)
{
    /*
        Function to performe mutation between 3 individuals.
                
                Parameters.
                                vector<vector<int>> individuals
                                float mutpb  -> mutation probability
    */
    
    int genMutating = 3, i, j, ind;
    std::uniform_int_distribution<> rand_num(0, individuals[0].size() - 1);
    
    for(i = 0; i < individuals.size(); i++)
    {
        for(j = 0; j < genMutating; j++)
        {
            if(dis(gen) < PARAMETERS["MUTATION_PROBABILITY"])
            {
                ind = rand_num(gen);
                individuals[i][ind] = 1 - individuals[i][ind];
            }
        }
    }
}


void reArrangeBits(std::vector<std::vector<int>> individuals)
{
    /*
        Function to performe mutation of all genes. We will rearrange
        the bits from  back to front. For example:
                {1 0 1 0} -> {0 1 0 1}
                
                Parameters.
                                vector<vector<int>> individuals
                                float mutpb  -> mutation probability
    */
    
    int i;
    
    for(i = 0; i < individuals.size(); i++)
    {
        // Reaarange from back to front
        std::reverse(individuals[i].begin(), individuals[i].end());
    }
}


/* ************************************************************************ */
