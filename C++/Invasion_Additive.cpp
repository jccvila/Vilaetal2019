// defaulte header files
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>      /* printf */
#include <math.h>       /*   */
#include <algorithm>
#include <cmath>
#include <cstdlib> // for exit()
#include <string>
#include <time.h>
#include <ctime>
#include <set>
#include <numeric>
#include <gsl/gsl_rng.h>


using namespace std;


/************************************************************
Default functions that deal with compiler version issues.
 ************************************************************/

template <typename T>
std::string to_string(T value)
{
    std::ostringstream os ;
    os << value ;
    return os.str() ;
}


template <class ForwardIterator, class T>
  void myiota (ForwardIterator first, ForwardIterator last, T val)
{
  while (first!=last) {
    *first = val;
    ++first;
    ++val;
  }
}

/************************************************************
Random number generator objects
 ************************************************************/
class NRrand

{
private:
  long s;
  const gsl_rng_type * T;
  gsl_rng * r;
  double u;
public:


void set_seed(long seed){
    s = seed;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set (r, s);
}

double d01()
{
  u = gsl_rng_uniform (r);
  return u;
}

long i0(long max)
{
    //long temp = (long(d01()*(max+1)));
    //cout << "random call = " << temp << "\n";
    return (long(d01()*(max+1)));
}

void clear_data(){
  gsl_rng_free (r);
}


}
;

/*
class NRrand
{
public:
    float d01(){return((float) rand() / (RAND_MAX));}
    int i0(int max){return((rand() % max));}
}
;
*/
/************************************************************
Simulation Object (basically main object)
 ************************************************************/
class Simulation
{
private:

    // list of core paramaters
    float mutation_rate;
    float selection_strength;
    int propagule_pressure;
    int task_number;
    int seed;

    // population and propaguel features to calculate over each run
    vector<int> propagule_fitness_mean;
    vector<int> propagule_fitness_variance;
    vector<int> propagule_richness;

    vector<float> invaded_mean_at_start;
    vector<float> invaded_variance_at_start;
    vector<int> invaded_richness_at_start;
    vector<float> invador_mean_at_start;
    vector<float> invador_variance_at_start;
    vector<int> invador_richness_at_start;

    vector<float> invaded_mean_at_end;
    vector<float> invaded_variance_at_end;
    vector<int> invaded_richness_at_end;
    vector<float> invador_mean_at_end;
    vector<float> invador_variance_at_end;
    vector<int> invador_richness_at_end;

    // calculated in each run

    vector<float> invasion_frequency;
    vector<int> invasion_fitness_invador;
    vector<int> invasion_fitness_invaded;

    // paramaters for home
    int home_popsize ;// home popsize
    int home_generations; // number of generations.

    // paramaters for evolution before invasion
    vector<int> invador_popsize ;// home popsize
    int temp_invador_popsize;
    vector<int> invador_generations; // number of generations.
    int temp_invador_generations;

    // paramatesr for post invasion
    int invaded_popsize ;// home popsize
    int invaded_generations; // number of generations.

    // directory paths
    string DirectoryPath;

    // home populations
    vector<int> home_population;  // starting population.
    vector<int> home_types;
    vector<int> home_type_abundances;
    vector<int> home_type_categories;

    // invading populations
    vector<int> invador_population;  // starting population.
    vector<int> invador_types;
    vector<int> invador_type_abundances;
    vector<int> invador_type_categories;

    // invadied populations
    vector<int> invaded_population;  // starting population.
    vector<int> invaded_types;
    vector<int> invaded_type_abundances;
    vector<int> invaded_type_categories;

    // used by dilution expansion step
    vector<int> new_population;
    vector<int> new_types;
    vector<int> new_type_categories;
    vector<int> new_type_abundances;

    //and the special ones for post invasion
    vector<bool> invasion_tracker;

    NRrand rng; // random number generator that is stored in header file.
    // file vectors

    // lists of more temporary filed objects used for saving and simulation. Swap between temporary and fixed objects to control simulation flow
    vector<int> population;
    vector<int> types;
    vector<int> type_abundances;
    vector<int> type_categories;

    // int recyler and type related objects

    // population level traits.
    int popsize;
    int min_cat;
    int max_cat;
    float invfreq;
    int invfitness_invador;
    int invfitness_invaded;


    // checkers
    bool check_max_cat ;
    bool check_min_cat ;
    // tracks stage of simulation. 0 is burnin, 1 is invador evolution, 2 is invaded evolution and 3 is the invasion.
    int check_stage;
    // updates paramatesr
    int update_params_counter;
    int param_run_number;
    float previous_runtime;

    // used types vector.
    vector <int> recycler;

    // stuff used for saving
    string temporary_DirectoryPath;
    vector<float> floats;
    vector<int> ints;
    // vector 0of ints is always a vector invasion invasion sums, vector of floats is a vector invasion freq.
    vector<vector<int> > vectors_of_ints_invaded;
    //invaded
    vector<vector<int> > vectors_of_ints_invador;
    //invaded
    vector<vector<float> > vectors_of_floats;
    // stuff used for iterating;
    // objects for iterations
    vector<int> popsizes;
    vector<int> generations;

    // indices used in by invasiona and dilution step
    vector<int> indices;
    vector<int> indices_a;
    vector<int> indices_b;


public:
    Simulation() {} // the most basic possible constructor

    /** Functions that make up generationstep **/

    float fitness_calculator(int fitness_category)
    {
        double fitness = 1 + (selection_strength * fitness_category);
        return fitness;
    }

    void reset_recycler()
    {
        //reemove start and end of type_abundances
        recycler.clear();
        while(type_abundances.front() == 0)
                { // i.e a type at lower end has gone extinct delete it, update min _type.
                type_abundances.erase(type_abundances.begin());
                type_categories.erase(type_categories.begin());
                transform(types.begin(),types.end(),types.begin(),bind2nd(std::plus<int>(),-1));
                }

        while(type_abundances.back() == 0){
                type_abundances.pop_back();
                type_categories.pop_back();
            }
        //reeset recycler
        for(unsigned int i = 0; i <type_abundances.size(); i++)
        {
            if(type_abundances[i] == 0)
            {
                recycler.push_back(i);
                type_categories[i] = -1;
            }
        }
    }
    // takes index born and mutate, also update types.
    void mutate(int index)
    {

        if (rng.d01() < mutation_rate) // mutation happens
        {
            // 50-50 change of up or down
            if (rng.d01() > 0.5)
            {
                    // if new individiaul is maximum update max cat. If they are minimum check min_cat as might have mutated away from minimum. Conservative test as multiple types may have minimum
                    if(population[index] == max_cat){max_cat++;}
                    if(population[index] == min_cat && type_abundances[types[index]] == 1){check_min_cat = true;}
                    population[index]++;
            }

            else
            {
                    // if new individiaul is maximum check max cat if minimum check min cat as might have mutated away from maximum. Conservative test as multiple types may have maximum
                    if(population[index] == min_cat){min_cat--;}
                    if(population[index] == max_cat && type_abundances[types[index]] == 1){check_max_cat = true;}
                    population[index]--;
            }

            // update original type. If types abundance is 0 add to recycler and set type category to -1
            type_abundances[types[index]]--;
            if(type_abundances[types[index]] == 0){
                            recycler.push_back(types[index]);
                            type_categories[types[index]] = -1;
                            }
            // update new type. If recycler is empty create new type else assign the new mutant the type at the end of the recycler.
            // also update categories by adding adding category to vector i fnew type or update category if it's a recycled type.

            if(recycler.empty()){
                types[index] = type_categories.size();
                type_categories.push_back(population[index]);
                type_abundances.push_back(1);
                }
            else{
                types[index] = recycler.front();
                type_categories[types[index]] = population[index];
                type_abundances[types[index]]++;
                recycler.erase(recycler.begin());
                }
        }

    }


    // takes population, and randomly picks an index.  That indvidual at that index is dead
    int death()
    {
        int index_a = rng.i0(popsize-1);
        return index_a;
    }

    // takes population and selection strength, picks an index, that individual is born if meets fitness criteria if not repeat.
    int birth()
    {


        float max_fitness = fitness_calculator(max_cat);
        int index_b;
        float fitness;
        float limit = rng.d01();

        // rng is random number simulator that pixs a random number, calculate it's fitness  and
        // then tests if fitness relative to max is < the randomly sampled number. This is a simple and quick method 
	// for implementing selection without having to resolve whole fitness distribution for every individual.

        do
        {
            index_b = rng.i0(popsize-1);
            fitness = fitness_calculator(population[index_b]);
        }
        while (fitness / max_fitness < limit);
        return index_b;
    }

    // performs n/2 rounds of bird death steps on population to count as singel generation
    // remember it's overlappping generations.
    void generationstep()
    {
        int index_a;
        int index_b;

        for( int a = 0; a < popsize/2; a++){
            do
            {
                index_a = death();
                index_b = birth();
            }
            while (index_a == index_b);
            // if the individual that died had the maximum or minimum category value and it's species is now extinct, it is worth checking if min cat or max cat
            // is correct. This is a conservative test as multiple species can have same category so only a subset of cases will min and max cat actually change.
            if(population[index_a] == min_cat && type_abundances[types[index_a]] == 1){// && types[dead_type-1] == 0){
                check_min_cat = true;
                }
            if(population[index_a] == max_cat && type_abundances[types[index_a]] == 1){//&& types[dead_type-1] == 0){
                check_max_cat = true;}
            // update abundances vectors and population and types vector.(no need to change type_category vector as only changes when new type is added or removed.
            type_abundances[types[index_a]]--;
            type_abundances[types[index_b]]++;
            types[index_a] = types[index_b];
            population[index_a] = population[index_b];
            if(check_stage == 3){invasion_tracker[index_a] = invasion_tracker[index_b];}
            // mutate function aslo determines if need to check min cat or max cat due to mutation.
            mutate(index_a);
            //update max and min category
            if(check_min_cat)
            {
                min_cat = *min_element(population.begin(), population.end());
                check_min_cat = false;
            }
            if(check_max_cat){
                max_cat = *max_element(population.begin(), population.end());
                check_max_cat = false;
            }


        }

    }

    /** Functions that save results to file**/

    // function that saves paramaters
    /** Functions that save results to file**/

    // function that saves paramaters
    void saveparams(string &FilePath)
    {
        // saving paramaters to file
        ofstream output; // initialise outf stream for paramaters
        const char * Path = FilePath.c_str();
        // convert const char
        output.open(Path,ios::binary);

        if( !output.is_open() )
        {
            cerr<<"\nERROR with file destination:["<<FilePath<<"]"<<endl;
            for(int i=0; i<5; i++)
            {
                cerr<<"********************************************************"<<endl;
            }
            cerr<<"******************PROGRAM TERMINATED********************\n";
            exit(1);
        }        // save paramaters
        output << "seed, " << seed << endl;
        output << "task_number, " << task_number/param_run_number << endl;
        output << "selection_strength, " << selection_strength << endl;
        output << "mutation_rate, " << mutation_rate << endl;
        output << "propagule_pressure, " << propagule_pressure << endl;
        output << "home_popsize, " << home_popsize << endl;
        output << "home_generations, " << home_generations << endl;
        output << "invador_popsize, " << temp_invador_popsize << endl;
        output << "invador_generations, " << temp_invador_generations << endl;
        output << "invaded_popsize, " << invaded_popsize << endl;
        output << "invaded_generations, " << invaded_generations << endl;
        output.close();
    }


    // function takes vecto of vectors of ints saves it into a .csv file
    void savevectorints_invaded(string &FilePath)
    {
        ofstream output; // initialise outf stream
        const char * Path = FilePath.c_str();
        // convert const char
        output.open(Path,ios::binary);
        if( !output.is_open() )
        {
            cerr<<"\nERROR with file destination:["<<FilePath<<"]"<<endl;
            for(int i=0; i<5; i++)
            {
                cerr<<"********************************************************"<<endl;
            }
            cerr<<"******************PROGRAM TERMINATED********************\n";
            exit(1);
        }
        int numberofrows= vectors_of_ints_invaded.size();
        for (int i = 0; i<= numberofrows- 1; i++)
        {
            ints = vectors_of_ints_invaded[i];
            int intsiterator = ints.size()-1;
            for (int i = 0; i <= intsiterator; i++)
            {
                output << ints[i];
                if (i < intsiterator)
                {
                    output << "," ;
                }

            }
            output<< endl;
            ints.clear();
        }
        vectors_of_ints_invaded.clear();
        output.close();
    }

    void savevectorints_invador(string &FilePath)
    {
        ofstream output; // initialise outf stream
        const char * Path = FilePath.c_str();
        // convert const char
        output.open(Path,ios::binary);
        if( !output.is_open() )
        {
            cerr<<"\nERROR with file destination:["<<FilePath<<"]"<<endl;
            for(int i=0; i<5; i++)
            {
                cerr<<"********************************************************"<<endl;
            }
            cerr<<"******************PROGRAM TERMINATED********************\n";
            exit(1);
        }
        int numberofrows= vectors_of_ints_invador.size();
        for (int i = 0; i<= numberofrows- 1; i++)
        {
            ints = vectors_of_ints_invador[i];
            int intsiterator = ints.size()-1;
            for (int i = 0; i <= intsiterator; i++)
            {
                output << ints[i];
                if (i < intsiterator)
                {
                    output << "," ;
                }

            }
            output<< endl;
            ints.clear();
        }
        vectors_of_ints_invador.clear();
        output.close();
    }

    void savevectorfloats(string &FilePath)
    {
        ofstream output; // initialise outf stream
        const char * Path = FilePath.c_str();
        // convert const char
        output.open(Path,ios::binary);
        if( !output.is_open() )
        {
            cerr<<"\nERROR with file destination:["<<FilePath<<"]"<<endl;
            for(int i=0; i<5; i++)
            {
                cerr<<"********************************************************"<<endl;
            }
            cerr<<"******************PROGRAM TERMINATED********************\n";
            exit(1);
        }
        int numberofrows= vectors_of_floats.size();
        for (int i = 0; i<= numberofrows- 1; i++)
        {
            floats = vectors_of_floats[i];
            int floatsiterator = floats.size()-1;
            for (int i = 0; i <= floatsiterator; i++)
            {
                output << floats[i];
                if (i < floatsiterator)
                {
                    output << "," ;
                }

            }
            output<< endl;
            floats.clear();
        }
        vectors_of_floats.clear();
        output.close();
    }

    void savesummarydata(string &FilePath)
    {
        // all vectors are the same size so generic iterator length.
        int vector_iterator =  param_run_number-1;
        ofstream output; // initialise outf stream
        const char * Path = FilePath.c_str();
        // convert const char
        output.open(Path,ios::binary);
        if( !output.is_open() )
        {
            cerr<<"\nERROR with file destination:["<<FilePath<<"]"<<endl;
            for(int i=0; i<5; i++)
            {
                cerr<<"********************************************************"<<endl;
            }
            cerr<<"******************PROGRAM TERMINATED********************\n";
            exit(1);
        }
        output << "Propagule_Fitness_Mean,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << propagule_fitness_mean[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Propagule_Fitness_Variance,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << propagule_fitness_variance[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Propagule_Richness,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << propagule_richness[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invader_Mean_Start,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invador_mean_at_start[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invader_Variance_Start,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invador_variance_at_start[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invader_Richness_Start,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invador_richness_at_start[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
            output << "Invader_Mean_End,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invador_mean_at_end[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invader_Variance_End,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invador_variance_at_end[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invader_Richness_End,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invador_richness_at_end[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
                output << "Invaded_Mean_Start,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invaded_mean_at_start[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invaded_Variance_Start,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invaded_variance_at_start[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invaded_Richness_Start,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invaded_richness_at_start[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invaded_Mean_End,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invaded_mean_at_end[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invaded_Variance_End,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invaded_variance_at_end[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output << "Invaded_Richness_End,";
        for (int i = 0; i <= vector_iterator; i++)
        {
            output << invaded_richness_at_end[i];
            if (i < vector_iterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        output.close();

        propagule_fitness_mean.clear();
        propagule_fitness_variance.clear();
        propagule_richness.clear();
        invaded_mean_at_end.clear();
        invaded_mean_at_start.clear();
        invaded_variance_at_end.clear();
        invaded_variance_at_start.clear();
        invaded_richness_at_end.clear();
        invaded_richness_at_start.clear();

        invador_mean_at_end.clear();
        invador_mean_at_start.clear();
        invador_variance_at_end.clear();
        invador_variance_at_start.clear();
        invador_richness_at_end.clear();
        invador_richness_at_start.clear();

    }
    /** Functions that calculate an output at intervals**/


    // calculator mean fitness category of population
    float fitness_category_mean_calculator()
    {
        float sum = 0;
        for(int i = 0; i < popsize; i++)
        {
            sum += population[i];
        }
        return sum/popsize;
    }

    // calculates fitness variance of population
    float fitness_category_variance_calculator()
    {
        float sum = 0;
        float sum_squared = 0;
        for(int i = 0; i < popsize; i++)
        {
            sum += population[i];
            sum_squared += pow(population[i],2);
        }
        return ((sum_squared) - (pow(sum,2)) / popsize) /popsize;
    }

    int richness_calculator()
    {
        int type_number = 0;
        for(unsigned int i = 0; i < type_abundances.size(); i++)
        {
            if(type_abundances[i] > 0){type_number++;}
        }
        return type_number;
    }
    void invador_frequency_calculator()
    {
        int invador_number = 0;
        invfreq = 0;
        invfitness_invaded = 0;
        invfitness_invador = 0;
        for (unsigned i=0; i < popsize; i++)
        {
            if (invasion_tracker[i] == true){
                invfitness_invador = invfitness_invador + population[i];
                invador_number++;
                }
            if(invasion_tracker[i] == false){
                invfitness_invaded = invfitness_invaded + population[i];
                }
        }
        invfreq = (float)invador_number / (float)popsize;
        // convert sum to mean correcting for proportion in population.
        if(invfreq != 0 & invfreq != 1){
            invfitness_invador = invfitness_invador/(invfreq * (float)popsize);
            invfitness_invaded = invfitness_invaded/((1-invfreq) * (float)popsize);
        }
        if(invfreq == 1){
            invfitness_invador = invfitness_invador/(float)popsize;
            invfitness_invaded = 0;
        }
        if(invfreq == 0){
            invfitness_invador = 0;
            invfitness_invaded = invfitness_invaded/(float)popsize;
        }
    //        cout << invador_frequency << endl;
    }

 /** Function performing manipulations**/

     void FisherYatesShuffle(vector<int> &indices){
        for (size_t k = 0; k < indices.size(); k++) {
            int r = k + rng.i0(indices.size()-k -1);
            swap(indices[k], indices[r]);
    //  std::cout << indices[r] << "\n";
            }
        }

    // currently worsk for both dilution and expansion of population (thoug diluations are only really needed_ does not alter fixed values only vectors.
    void dilutionexpansionstep(int newpopsize,bool check_invador)
    {
        // optain shuffled vector of indices of oldpopulation
        int oldpopsize = home_population.size();
        indices.resize(oldpopsize);
        int new_index; // for x`
        myiota(indices.begin(),indices.end(),0);
        FisherYatesShuffle(indices);
        // type_abundances are 0 and type categories are stored.
        new_type_abundances.assign(home_type_abundances.size(),0);
        new_type_categories = home_type_categories;
        // now create new objects
        if (newpopsize <= oldpopsize){
            for (int i = 0; i < newpopsize;i++){
                new_population.push_back(home_population[indices[i]]);
                new_types.push_back(home_types[indices[i]]);
                new_type_abundances[new_types[i]]++;
                if (new_type_abundances[new_types[i]] == 0){printf("type_abundance_error 2");}

            }
        }
        else
        {
            for (int i = 0; i < oldpopsize;i++)
            {
                new_population.push_back(home_population[indices[i]]);
                new_types.push_back(home_types[indices[i]]);
                new_type_abundances[new_types[i]]++;
                if (new_type_abundances[new_types[i]] == 0){printf("type_abundance_error 2");}
            }
            for (int i = oldpopsize; i < newpopsize;i++){
                new_index = rng.i0(oldpopsize - 1);
                new_population.push_back(home_population[new_index]);
                new_types.push_back(home_types[new_index]);
                new_type_abundances[new_types[i]]++;
                if (new_type_abundances[new_types[i]] == 0){printf("type_abundance_error 2");}

            }

        }

        // update global objects
        if(check_invador==true)
        {
            invador_population.clear();
            invador_types.clear();
            invador_type_abundances.clear();
            invador_type_categories.clear();
            invador_population = new_population;
            invador_types = new_types;
            invador_type_abundances = new_type_abundances;
            invador_type_categories = new_type_categories;
        }
        if(check_invador==false)
        {
            invaded_population.clear();
            invaded_types.clear();
            invaded_type_abundances.clear();
            invaded_type_categories.clear();
            invaded_population = new_population;
            invaded_types = new_types;
            invaded_type_abundances = new_type_abundances;
            invaded_type_categories = new_type_categories;

        }
        // clearer
        new_population.clear();
        new_types.clear();
        new_type_categories.clear();
        new_type_abundances.clear();
        indices.clear();
    }

    void invasion()
    {
        // optain shuffled vector of indices of oldpopulation
        int invador_pop_size  = invador_population.size();
        int invaded_pop_size  = invaded_population.size();
        int temp_propagule_fitness = 0;
        int temp_propagule_fitness_squared = 0;
        int index = 0;
        indices_a.resize(invador_pop_size);
        indices_b.resize(invaded_pop_size);
        myiota(indices_a.begin(),indices_a.end(),0);
        myiota(indices_b.begin(),indices_b.end(),0);
        FisherYatesShuffle(indices_a);
        FisherYatesShuffle(indices_b);
        vector<int> old_types_list;
        vector<int> new_types_list;


        // if for some weird reason propagule pressure is larger than population it originates from
        if(propagule_pressure > invador_pop_size)
        {
            for(int i = 0; i < (propagule_pressure - invador_pop_size);i++)
            {
                indices_a.push_back(rng.i0(invador_pop_size - 1));
                printf("prop pressure bigger than popsize");
            }
        }
        // and same below
        if(propagule_pressure > invaded_pop_size)
        {
            for(int i = 0; i < (propagule_pressure - invaded_pop_size);i++)
            {
                indices_b.push_back(rng.i0(invaded_pop_size - 1));
            }
        }

        for (int i = 0; i < propagule_pressure;i++)
        {
            //individual replae
            invaded_type_abundances[invaded_types[indices_b[i]]]--;
            if( invaded_type_abundances[invaded_types[indices_b[i]]] < 0){printf("invasion error 1");}
            invaded_population[indices_b[i]] = invador_population[indices_a[i]];
            temp_propagule_fitness = temp_propagule_fitness + invador_population[indices_a[i]];
            temp_propagule_fitness_squared = temp_propagule_fitness_squared + pow(invador_population[indices_a[i]],2);

            invasion_tracker[indices_b[i]] = true;

            // old types is'nt empty and invador type is found in older types list.
            if(!old_types_list.empty() && std::find(old_types_list.begin(), old_types_list.end(), invador_types[indices_a[i]]) != old_types_list.end())
            {
                index = std::find(old_types_list.begin(), old_types_list.end(), invador_types[indices_a[i]]) - old_types_list.begin();
                invaded_types[indices_b[i]] = new_types_list[index];
                invaded_type_abundances[new_types_list[index]]++;
            }
            else
            {
                // trakc old and new types and then put them into list
                old_types_list.push_back(invador_types[indices_a[i]]);
                new_types_list.push_back(invaded_type_categories.size());
                invaded_types[indices_b[i]] = invaded_type_categories.size();
                invaded_type_abundances.push_back(1);
                invaded_type_categories.push_back(invaded_population[indices_b[i]]);
            }
        }
        propagule_richness.push_back(old_types_list.size());
        propagule_fitness_mean.push_back((float) temp_propagule_fitness /(float)propagule_pressure);
        propagule_fitness_variance.push_back(((float) temp_propagule_fitness_squared - (float) pow(temp_propagule_fitness,2) /(float)propagule_pressure ) /float(propagule_pressure) );

        // reset types.
        // subtract min_type from every type.
        indices_a.clear();
        indices_b.clear();
        old_types_list.clear();
        new_types_list.clear();


    }
    /** Function that runs a neutral theory on "population" and saves results at intervals**/

    // runs simulation for seed, generation, strenght, mutation, rate population, pop interval, can specify directory path in which to ave.
    void runsimulation(int generations)
    {
        reset_recycler();
        // if(type_abundances.size() != type_categories.size()){printf("pre sim error 1");}
        // if(types.size() != population.size()){printf("pre sim error 2");}
        // if(population.size() != (unsigned)std::accumulate(type_abundances.begin(), type_abundances.end(), 0)){printf("pre sim error 3");}
        // for(unsigned int i = 0; i < types.size(); i++)if(types[i] < 0){printf("pre sim error 4");}


        // set values
        popsize = population.size();
        min_cat = *min_element(type_categories.begin(), type_categories.end());
        max_cat = *max_element(type_categories.begin(), type_categories.end());
        invfreq = 0.5; // to not quit loop on othe runs

        if (check_stage == 1)
        {
            invador_mean_at_start.push_back(fitness_category_mean_calculator());
            invador_variance_at_start.push_back(fitness_category_variance_calculator());
            invador_richness_at_start.push_back(richness_calculator());
        }
        if (check_stage == 2)
        {
            invaded_mean_at_start.push_back(fitness_category_mean_calculator());
            invaded_variance_at_start.push_back(fitness_category_variance_calculator());
            invaded_richness_at_start.push_back(richness_calculator());

        }
        // if it's an invasion fill in details

        if(check_stage == 3)
        {
            // by default
            invfitness_invaded = 0;
            invfitness_invador = 0;
            invfreq = 0;
            invador_frequency_calculator();
            invasion_frequency.push_back(invfreq);
            invasion_fitness_invaded.push_back(invfitness_invaded);
            invasion_fitness_invador.push_back(invfitness_invador);
        }

        // runs simulation update outputs.
        for( int a = 1; a <= generations; a = a + 1 )
        {
            generationstep();
            reset_recycler();
            if(check_stage == 3)
            {
                // by default
                invador_frequency_calculator();
                invasion_frequency.push_back(invfreq);
                invasion_fitness_invaded.push_back(invfitness_invaded);
                invasion_fitness_invador.push_back(invfitness_invador);
            }
            if(invfreq == 0 || invfreq == 1)break;
        }

        //output same depends on which run i'm doing. '
        if (check_stage == 1)
        {
            invador_mean_at_end.push_back(fitness_category_mean_calculator());
            invador_variance_at_end.push_back(fitness_category_variance_calculator());
            invador_richness_at_end.push_back(richness_calculator());
        }
        if (check_stage == 2)
        {
            invaded_mean_at_end.push_back(fitness_category_mean_calculator());
            invaded_variance_at_end.push_back(fitness_category_variance_calculator());
            invaded_richness_at_end.push_back(richness_calculator());

        }
        if (check_stage == 3)
        {
            vectors_of_floats.push_back(invasion_frequency);
            vectors_of_ints_invaded.push_back(invasion_fitness_invaded);
            vectors_of_ints_invador.push_back(invasion_fitness_invador);
            invasion_frequency.clear();
            invasion_fitness_invaded.clear();
            invasion_fitness_invador.clear();

        }

       // if(type_abundances.size() != type_categories.size()){printf("post sim error 1");}
       // if(types.size() != population.size()){printf("post sim error 2");}
       // if(population.size() != (unsigned)std::accumulate(type_abundances.begin(), type_abundances.end(), 0)){printf("post sim error 3");}
       // for(unsigned int i = 0; i < types.size(); i++)if(types[i] < 0){printf("post sim error 4");}
    }
    /** Function that sets initail values for private objects**/

    // setter, sets all the private variables.
    void setup(int internal_seed_number, float internal_selection_strength, float internal_mutation_rate,
        int internal_propagule_pressure, int internal_home_popsize)
    {

        //set paramaters that are parsed
        seed = internal_seed_number;
        selection_strength =internal_selection_strength;
        mutation_rate = internal_mutation_rate;
        propagule_pressure = internal_propagule_pressure;

         // internally defined and varied
        task_number = 0;

        // construct vectors for internally varied paramaters
        popsizes.clear();
        generations.clear();

        popsizes.push_back(500);
        popsizes.push_back(1000);
        popsizes.push_back(2000);
        popsizes.push_back(3000);
        popsizes.push_back(6000);
        popsizes.push_back(10000);
        popsizes.push_back(20000);
        popsizes.push_back(30000);
        popsizes.push_back(60000);
        popsizes.push_back(100000);

        generations.push_back(50);
        generations.push_back(100);
        generations.push_back(200);
        generations.push_back(350);
        generations.push_back(600);
        generations.push_back(1050);
        generations.push_back(1850);
        generations.push_back(3350);
        generations.push_back(6000);
        generations.push_back(10000);

        // now construct vectors from them to get all combinations
        for (std::size_t i = 0 ; i < popsizes.size() ; i ++)
        {
        for (std::size_t j = 0 ; j < generations.size() ; j ++)
        {
        invador_popsize.push_back(popsizes[i]);  // popsizes
        invador_generations.push_back(generations[j]); // number of generations.
        }
        }


        // home
        home_popsize = internal_home_popsize ;
        home_generations = home_popsize*4; // run for popsize number of generations for burnin

        // invaded
        invaded_popsize = home_popsize;
        invaded_generations = 100000000; // number of generations.

        // assign default values to populations
        home_population.assign(home_popsize,100);
        home_types.assign(home_popsize,0);
        home_type_abundances.assign(1,home_popsize);
        home_type_categories.assign(1,100);


        check_max_cat = false ;// default
        check_min_cat = false;// default
        check_stage = 0;
        // set seed
        rng.set_seed(seed);
        // paths
        DirectoryPath = "" ;

        // set up temporary objects
        population = home_population;
        types = home_types;
        type_abundances = home_type_abundances;
        type_categories = home_type_categories;

        // tracked values
        popsize = population.size();
        max_cat = 100;
        min_cat = 100;
        update_params_counter = 0;
        param_run_number = 100;
        if(home_popsize == 100000){param_run_number = 10;}
        if(home_popsize == 50000){param_run_number = 20;}
        if(home_popsize == 10000){param_run_number = 100;}
        if(home_popsize == 5000){param_run_number = 200;}
        if(home_popsize == 1000){param_run_number = 500;}


        previous_runtime = 0;
    }

    /** Function thats called to piece verything together**/

      void run(int internal_seed,float internal_selection_strength, float internal_mutation_rate, int internal_propagule_pressure,
                int internal_home_popsize, int wall_time)
    {

        // set up times
        time_t start , end;
        time(&start);
        bool timeout = false;
        setup(internal_seed, internal_selection_strength, internal_mutation_rate,
                internal_propagule_pressure, internal_home_popsize);
        /** burn in**/


        // run burn in
        cout << "Starting Simulation " << seed << endl;
        cout << "Propagule Pressure " << propagule_pressure << endl;
        cout << "Mutation Rate " << mutation_rate << endl;
        cout << "Selection Strengt " << selection_strength<< endl;
        cout << "Population Size " << home_popsize << endl;
        cout << "Burnin Begin" << endl;
        runsimulation(home_generations);

        // save tracked values for resetting later on
        //clear
        home_population.clear();
        home_types.clear();
        home_type_abundances.clear();
        home_type_categories.clear();
        // setup invasion and invaded population from this result

        home_population = population;
        cout << fitness_category_mean_calculator();
        home_types = types;
        home_type_abundances = type_abundances;
        home_type_categories =type_categories;
        cout << "Burnin End" << endl;

        // clear tempoary objects so that we're starting from correct state
        population.clear();
        types.clear();
        type_abundances.clear();
        type_categories.clear();
        // now for the loop

        do
        {
            /** taks allocation part **/
            task_number++;
            rng.clear_data();
            rng.set_seed((1000 * task_number) + seed);
            if(task_number % param_run_number == 1){
                temp_invador_generations = invador_generations[update_params_counter];
                temp_invador_popsize = invador_popsize[update_params_counter];
                update_params_counter++;
                }

            /** set up populations **/

            cout << "Starting task "  << (task_number-1)/param_run_number + 1  <<   " part "   << (task_number-1)%param_run_number  +1 << endl;
            // reset values from previous run ( for invador evolution)
            // create_new_populations
            dilutionexpansionstep(temp_invador_popsize,true); // set up invador (remember it varies by run
            dilutionexpansionstep(invaded_popsize,false); // set up invaded (doesnt vary by run.

            // error checer section 1
            if(home_type_abundances.size() != home_type_categories.size()){printf("home population error 1\n");}
            if(home_types.size() != home_population.size()){printf("home population error 2\n");}
            if(home_population.size() !=  (unsigned)std::accumulate(home_type_abundances.begin(), home_type_abundances.end(), 0)){printf("home population error 3\n");}
            for(unsigned int i = 0; i < home_types.size(); i++) if(home_types[i] < 0){printf("home population error 4\n");}

            if(invaded_type_abundances.size() != invaded_type_categories.size()){printf("invaded population error 1\n");}
            if(invaded_types.size() != invaded_population.size()){printf("invaded population error 2\n");}
            if(invaded_population.size() != (unsigned)std::accumulate(invaded_type_abundances.begin(), invaded_type_abundances.end(), 0)){printf("invaded population error 3\n");}
            for(unsigned int i = 0; i < invaded_types.size(); i++)if(invaded_types[i] < 0){printf("invaded population error 4\n");}

            if(invador_type_abundances.size() != invador_type_categories.size()){printf("invador population error 1\n");}
            if(invador_types.size() != invador_population.size()){printf("invador population error 2\n");}
            if(invador_population.size() != (unsigned)std::accumulate(invador_type_abundances.begin(), invador_type_abundances.end(), 0)){printf("invador population error 3\n");}
            for(unsigned int i = 0; i < invador_types.size(); i++)if(invador_types[i] < 0){printf("invador population error 4\n");}

            /** invader evolution **/

            // convert invader vector to default vector and clear invader
            population = invador_population;
            types = invador_types;
            type_abundances = invador_type_abundances;
            type_categories = invador_type_categories;
            temporary_DirectoryPath = DirectoryPath;
            invador_population.clear();
            invador_types.clear();
            invador_type_abundances.clear();
            invador_type_categories.clear();
            // run invador evolution
            printf("Invador Evolution Begin\n");
            check_stage = 1;
            runsimulation(temp_invador_generations);
            printf("Invador Evolution End\n");
            // convert default vector to invader vector and clear default
            invador_population = population;
            invador_types = types;
            invador_type_abundances = type_abundances;
            invador_type_categories = type_categories;
            population.clear();
            types.clear();
            type_abundances.clear();
            type_categories.clear();
            //Run Invaded Evolution

            /** invaded evolution**/

            population = invaded_population;
            types = invaded_types;
            type_abundances = invaded_type_abundances;
            type_categories = invaded_type_categories;
            temporary_DirectoryPath = DirectoryPath;
            invaded_population.clear();
            invaded_types.clear();
            invaded_type_abundances.clear();
            invaded_type_categories.clear();
            // run invador evolution
            printf("Invaded Evolution Begin\n");
            check_stage = 2;
            runsimulation(temp_invador_generations);
            printf("Invaded Evolution End\n");
            // convert default vector to invader vector and clear default
            invaded_population = population;
            invaded_types = types;
            invaded_type_abundances = type_abundances;
            invaded_type_categories = type_categories;
            population.clear();
            types.clear();
            type_abundances.clear();
            type_categories.clear();
            // clear temporary objects

            /** Invasion **/

            // reset for back to burn in levels
            invasion_tracker.assign(invaded_popsize,false); // by default no one is an invador.
            check_stage = 3;
            // error checker 2 invasion
            if(invador_type_abundances.size() != invador_type_categories.size()){printf("pre invasion error 1\n");}
            if(invador_types.size() != invador_population.size()){printf("pre invasion error 2\n");}
            if(invador_population.size() != (unsigned)std::accumulate(invador_type_abundances.begin(), invador_type_abundances.end(), 0)){printf("pre invasion error 3\n");}
            for(unsigned int i = 0; i < invaded_types.size(); i++)if(invaded_types[i] < 0){printf("pre invasion error 4\n");}
            invasion();
            if(invaded_type_abundances.size() != invaded_type_categories.size()){printf("post invasion error 1\n");}
            if(invaded_types.size() != invaded_population.size()){printf("post invasion error 2\n");}
            if(invaded_population.size() != (unsigned)std::accumulate(invaded_type_abundances.begin(), invaded_type_abundances.end(), 0)){printf("post invasion error 3\n");}
            for(unsigned int i = 0; i < invaded_types.size(); i++)if(invaded_types[i] < 0){printf("post invasion error 4\n");}

            // now reset
            population = invaded_population;
            types = invaded_types;
            type_abundances = invaded_type_abundances;
            type_categories = invaded_type_categories;
            temporary_DirectoryPath = DirectoryPath;

            // run invador evolution
            printf("Invasion Begin\n");
            runsimulation(invaded_generations);
            printf("Invasion End\n");

            // fixed objects
            invador_population.clear();
            invador_types.clear();
            invador_type_abundances.clear();
            invador_type_categories.clear();

            // simulation used objects
            population.clear();
            types.clear();
            type_categories.clear();
            type_abundances.clear();
            invasion_tracker.clear();
            if(task_number % param_run_number == 0){
                string paramaters_path = DirectoryPath + "Paramaters_" + to_string(seed) + "_" + to_string((task_number/param_run_number)) +".csv";
                string invasion_freq_path =  temporary_DirectoryPath + "Invasion_Frequency_" + to_string(seed) + "_" + to_string((task_number/param_run_number)) +".csv";
                string invasion_fitness_invador_path = temporary_DirectoryPath + "Invasion_Fitness_Invader_" + to_string(seed) + "_" + to_string((task_number/param_run_number)) +".csv";
                string invasion_fitness_invaded_path = temporary_DirectoryPath + "Invasion_Fitness_Invaded_" + to_string(seed) + "_" + to_string((task_number/param_run_number)) +".csv";
                string summary_data_path = temporary_DirectoryPath + "Summary_Data_" + to_string(seed) + "_" + to_string((task_number/param_run_number)) +".csv";
                saveparams(paramaters_path);
                savevectorfloats(invasion_freq_path);
                savevectorints_invaded(invasion_fitness_invaded_path);
                savevectorints_invador(invasion_fitness_invador_path);
                savesummarydata(summary_data_path);
            }

            cout << "Ending task " << task_number/param_run_number +1 <<  " part "  << task_number%param_run_number << endl;
            time(&end);
            cout << "Task took " <<  difftime(end,start) - previous_runtime << "second" << "\n"<< endl;
            previous_runtime = difftime(end,start);
//            if (difftime(end,start) >= (wall_time)){timeout = true;}
            if ((unsigned)task_number >= (invador_generations.size() * param_run_number)){
                cout << "simulation_complete" << endl;
                timeout = true;}
            if (difftime(end,start) >= wall_time){
                cout << "wall time hit" << endl;
                timeout = true;}
        }
        while (!timeout);
        cout << "Ending simulation " << seed << " at " << task_number<< " which took " <<  difftime(end,start) << endl;            // populations
    }

};

/************************************************************
 MAIN
 ************************************************************/

int charconvertor(char charin)
{
    switch (charin)
    {
    case '0':
        return (0);
    case '1':
        return (1);
    case '2':
        return (2);
    case '3':
        return (3);
    case '4':
        return (4);
    case '5':
        return (5);
    case '6':
        return (6);
    case '7':
        return (7);
    case '8':
        return (8);
    case '9':
        return (9);
    default:
        return (-1);
    }
}

int jobconvertor(char* argin)
{
    int maxind = 0;
    while (charconvertor(argin[maxind]) != -1)
    {
        maxind ++;
    }
    int jobtoret = 0;
    int pow10 = 1;
    for (int i = maxind-1 ; i >=0 ; i --)
    {
        jobtoret += (pow10*charconvertor(argin[i]));
        pow10 = pow10*10;
    }
    return jobtoret;
}




int main(int argc, char* argv[])
{

    // first get the command line arg data in

    int seed_number = -1;// the task number
    float wall_time = -1;// the max time (seconds)
//
////    // read command line data
    seed_number = jobconvertor(argv[1]);
    wall_time = jobconvertor(argv[2])  * 3600; //(convert hours to minute) ;
    vector<int> p_vect; // propagule pressure
    vector<int> jm_vect;
    vector<float> nu_vect; // mutation rates
    vector<float> s_vect; // selection strength

    p_vect.clear();
    jm_vect.clear();
    nu_vect.clear();
    s_vect.clear();


    // uninteresting paramater
    jm_vect.push_back(1000);
    jm_vect.push_back(5000);
    jm_vect.push_back(10000);
    jm_vect.push_back(50000);
    jm_vect.push_back(100000);

    nu_vect.push_back(0);
    nu_vect.push_back(0.00001);
    nu_vect.push_back(0.0001);

    s_vect.push_back(0);
    s_vect.push_back(0.001);
    s_vect.push_back(0.01);

    p_vect.push_back(1);
    p_vect.push_back(2);
    p_vect.push_back(3);
    p_vect.push_back(5);
    p_vect.push_back(10);
    p_vect.push_back(20);
    p_vect.push_back(30);
    p_vect.push_back(50);
    p_vect.push_back(100);
    p_vect.push_back(200);



    vector<int> p_list;
    vector<int> jm_list;
    vector<float> nu_list;
    vector<float> s_list;

    // construct paramater combinations
    for (unsigned int i = 0 ; i < jm_vect.size() ; i ++)
    {
    for (unsigned int j = 0 ; j < nu_vect.size() ; j ++)
    {
    for (unsigned int k = 0 ; k < s_vect.size() ; k ++)
    {
    for (unsigned int l = 0 ; l < p_vect.size() ; l ++)
    {
    p_list.push_back(p_vect[l]);
    jm_list.push_back(jm_vect[i]);
    nu_list.push_back(nu_vect[j]);
    s_list.push_back(s_vect[k]);
    }
    }
    }
    }

    cout << seed_number << " = seed \n";
    cout << wall_time << " = max time \n";
    Simulation simulation;
    simulation.run(seed_number,s_list[seed_number-1],nu_list[seed_number-1],p_list[seed_number-1],jm_list[seed_number-1],wall_time);

    return 0;

}


