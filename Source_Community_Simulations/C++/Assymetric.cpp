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

// custome header files
//#include "rng2.h" // self defined header to have rng
// definitions for rng
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 5277
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-8
#define RNMX (1.0-EPS)

using namespace std;


/************************************************************
Compiler fucked up functions
 ************************************************************/

// long maxval(vector<long>& vect){
//    long maximum;
//    maximum = vect[0];
//    for (unsigned int i = 1; i <- vect.size(); i++) if (vect[i] < maximum){ maximum = vect[i];}
//        return(maximum);
//    }
////
//long minval(vector<long>& vect){
//    long minimum = vect[0];
//    for (unsigned int i = 1; i < vect.size(); i++) if (vect[i] > minimum){minimum = vect[i];}
//        return(minimum);
//    }

template <typename T>
std::string to_string(T value)
{
    std::ostringstream os ;
    os << value ;
    return os.str() ;
}


/************************************************************
Random number generator objects
 ************************************************************/
class NRrand
{

private:
    long idum;
    int j;
    long k;
    long idum2;
    long iy;
    long iv[NTAB];
    double temp;
    bool seeded;

    double lastresult;
    bool normflag;


public:

    NRrand()
    {
        seeded = false;
        normflag = true;

    }

    void set_seed(long seed)
    {
        if (!seeded)
        {
            idum2 = 123456789;
            iy = 0;
            idum = seed;
            if (idum < 1) idum=1;
            //Be sure to prevent idum = 0.
            idum2=(idum);
            for (j=NTAB+7; j>=0; j--)
            {
                //Load the shuffle table (after 8 warm-ups).
                k=(idum)/IQ1;
                idum=IA1*(idum-k*IQ1)-k*IR1;
                if (idum < 0) idum += IM1;
                if (j < NTAB) iv[j] = idum;
            }
            iy=iv[0];
            seeded = true;
        }
    }

    double d01()
    {
        k=(idum)/IQ1;
        //Start here when not initializing.
        idum=IA1*(idum-k*IQ1)-k*IR1;
        //Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method.
        if (idum < 0) idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        //Compute idum2=(IA2*idum) % IM2 likewise.
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        //Will be in the range 0..NTAB-1.
        iy=iv[j]-idum2;
        //Here idum is shuffled, idum and idum2 are combined to generate output.
        iv[j] = idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX)
        {
            //cout << "random call = " << "RNMAX" << "\n";
            return RNMX; //Because users don't expect endpoint values.
        }
        else
        {

            return temp;
        }
    }

    long i0(long max)
    {
        //long temp = (long(d01()*(max+1)));
        //cout << "random call = " << temp << "\n";
        return (long(d01()*(max+1)));
    }

};

/*
class NRrand
{
public:
    double d01(){return((double) rand() / (RAND_MAX));}
    long i0(long max){return((rand() % max));}
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
    double mutation_rate;
    double selection_strength;
    long seed;
    long task_number;

    // paramaters for home
    long home_popsize ;// home popsize
    long home_generations; // number of generations.
    long home_interval_mean; // interval at which mean is measured
    long home_interval_variance;// interval at which variance is measured
    long home_interval_abundances; // interval at which abundance is measured

    // list of specific objects (i.e objects meant to appy to one of the three simulations
    // directory paths
    string DirectoryPath;

    // populations
    vector<long> home_population;  // starting population.
    vector<long> home_types;
    vector<long> home_type_abundances;
    vector<long> home_type_categories;
    NRrand rng; // random number generator that is stored in header file.
    // file vectors

    // lists of more temporary filed objects used for saving and simulation. Swap between temporary and fixed objects to control simulation flow
    vector<long> population;
    vector<long> types;
    vector<long> type_abundances;
    vector<long> type_categories;

    // long recyler and type related objects
    vector <long> recycler;
    long extinct_type_counter;

    // population level traits.

    long popsize;
    long min_cat;
    long max_cat;
    bool check_max_cat ;
    bool check_min_cat ;



    // stuff used for saving
    string temporary_DirectoryPath;
    vector<double> doubles;
    vector<long> longs;
    vector<vector<long> > vectors_of_longs;
    vector<vector<long> > vectors_of_doubles;

    vector<double> fitness_means; // vector of mean fitness(not categories
    vector<double> fitness_variance; // vector of fitness variances (not categories_)
    vector<double> fitness_category_means; // vector of mean fitness(not categories
    vector<double> fitness_category_variance; // vector of fitness variances (not categories_)
    vector<long>   extinct_counter; // vector of number of non tracked types
    vector<vector<long> >  vector_of_type_categories; // vector of type categories
    vector<vector<long> >  vector_of_type_abundances; // vector of type abundancs.

//    long counter = 1;
    long counter;
    long min_type;



public:
    Simulation() {} // the most basic possible constructor

    /** Functions that make up generationstep **/

    double fitness_calculator(long fitness_category)
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
                transform(types.begin(),types.end(),types.begin(),bind2nd(std::plus<long>(),-1));
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
    // takes index born and mutate, also update types.
    void mutate(int index)
    {

        if (rng.d01() < mutation_rate) // mutation happens
        {

            // 50-50 change of up or down
            if (rng.d01() > 0.9)
            {
                    // do stuff to population.
                    if(population[index] == max_cat){max_cat++;}
                    if(population[index] == min_cat && type_abundances[types[index]] == 1){check_min_cat = true;}
                    population[index]++;
            }

            else
            {
                    // do opposite to population
                    if(population[index] == min_cat){min_cat--;}
                    if(population[index] == max_cat && type_abundances[types[index]] == 1){check_max_cat = true;}
                    population[index]--;
            }



            // update types
            type_abundances[types[index]]--;
            if(type_abundances[types[index]] == 0){
                            recycler.push_back(types[index]);
                            extinct_type_counter++;
                            type_categories[types[index]] = -1;
                            }
            if(recycler.empty()){
                types[index] = type_categories.size();
                type_categories.push_back(population[index]);
                type_abundances.push_back(1);

                }
            else{
                types[index] = recycler[0];
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


        double max_fitness = fitness_calculator(max_cat);
        int index_b;
        double fitness;
        double limit = rng.d01();

        // rng is random number simulator this pixs a random number, calculate it's fitness (absolute value again) and
        // then tests of relative fitness to max is < randomly sampled number. This is a simple and quick to do selection
        // without having to resolve whole fitness distribution for every individual.
//        int loop_counter = 0;
        do
        {
            index_b = rng.i0(popsize-1);
            fitness = fitness_calculator(population[index_b]);
//            if (loop_counter = 500){
//                max_cat = *max_element(population.begin(), population.end());
//                max_fitness = fitness_calculator(max_cat);
//                }
//            loop_counter++;
//               cout << limit << endl;
            // std::cout << "\n"<< fitness << " divided by " << max_fitness << " equals " << fitness / max_fitness;
        }
        while (fitness / max_fitness < limit);
        return index_b;
    }

    // performs n/2 rounds of bird death steps on population to count as singel generation
    // remember it's overlappping generations.
    void generationstep()
    {
        long index_a;
        long index_b;
        counter =0;
        for( int a = 1; a < popsize/2; a = a + 1 )
            {
            do
            {
                index_a = death();
                index_b = birth();

            }
            while (index_a == index_b);
            // ok now do everything first check if max or min cat type might be extinct (not 100% of time but reduces risk of error.
            if(population[index_a] == min_cat && type_abundances[types[index_a]] == 1){// && types[dead_type-1] == 0){
                check_min_cat = true;
                }
            if(population[index_a] == max_cat && type_abundances[types[index_a]] == 1){//&& types[dead_type-1] == 0){
                check_max_cat = true;}
            // update abundances vector
            type_abundances[types[index_a]]--;
            type_abundances[types[index_b]]++;
            // no need to update max and min_cat except if they a) the individual dying has that value and b) their type is extinct.
            types[index_a] = types[index_b];
            population[index_a] = population[index_b];


            // never a need to update max_type ever used or to change the type to categries vector as no new types are originating.

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
    void saveparams()
    {
        // saving paramaters to file
        ofstream output; // initialise outf stream for paramaters
        string FilePath = DirectoryPath + "Paramaters_" + to_string(seed) + "_" + to_string(task_number) +".csv"; // define path
        const char * Path = FilePath.c_str();
        // convert const char
        output.open(Path,ios::binary);
        // save paramaters
        output << "seed, " << seed << endl;
        output << "task_number, " << task_number << endl;
        output << "selection_strength, " << selection_strength << endl;
        output << "mutation_rate, " << mutation_rate << endl;
        output << "home_popsize, " << home_popsize << endl;
        output << "home_generations, " << home_generations << endl;
        output << "home_interval_mean, " << home_interval_mean << endl;
        output << "home_interval_variance, " << home_interval_variance << endl;
        output << "home_interval_abundances, " << home_interval_abundances << endl;

        output.close();
    }

    // function takessaves current vector of doubles to location
    void savedoubles(string &FilePath)
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
        int doublesiterator = doubles.size() -1;
        for (int i = 0; i <= doublesiterator; i++)
        {
            output << doubles[i] ;
            if (i < doublesiterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        doubles.clear();
        output.close();
    }

        // function takessaves current vector of doubles to location
    void savelongs(string &FilePath)
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
        int longsiterator = longs.size() -1;
        for (int i = 0; i <= longsiterator; i++)
        {
            output << longs[i] ;
            if (i < longsiterator)
            {
                output << "," ;
            }
        }
        output<< endl;
        longs.clear();
        output.close();
    }

    // function takes vecto of vectors of longs saves it into a .csv file
    void savevectorlongs(string &FilePath)
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
        int numberofrows= vectors_of_longs.size();
        for (int i = 0; i<= numberofrows- 1; i++)
        {
            longs = vectors_of_longs[i];
            int longsiterator = longs.size()-1;
            for (int i = 0; i <= longsiterator; i++)
            {
                output << longs[i];
                if (i < longsiterator)
                {
                    output << "," ;
                }

            }
            output<< endl;
            longs.clear();
        }
        vectors_of_longs.clear();
        output.close();
    }

    void savevectordoubles(string &FilePath)
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
        int numberofrows= vectors_of_doubles.size();
        for (int i = 0; i<= numberofrows- 1; i++)
        {
            longs = vectors_of_doubles[i];
            int doublesiterator = doubles.size()-1;
            for (int i = 0; i <= doublesiterator; i++)
            {
                output << doubles[i];
                if (i < doublesiterator)
                {
                    output << "," ;
                }

            }
            output<< endl;
            doubles.clear();
        }
        vectors_of_doubles.clear();
        output.close();
    }
    /** Functions that calculate an output at intervals**/

    // calculator mean fitness of population
    double fitness_mean_calculator()
    {
        double sum = 0;
        double fitness;
        for(int i = 0; i < popsize; i++)
        {
            if (population[i] > 0)
            {
                fitness = fitness_calculator(population[i]);
            }
            else
            {
                fitness = fitness_calculator(population[i] * -1) ;
            }
            sum += fitness;
        }
        return sum/popsize;
    }

    // calculator mean fitness category of population
    double fitness_category_mean_calculator()
    {
        double sum = 0;
        for(int i = 0; i < popsize; i++)
        {
            sum += population[i];
        }
        return sum/popsize;
    }

    // calculates fitness variance of population
    double fitness_variance_calculator()
    {
        double sum = 0;
        double sum_squared = 0;
        double fitness;
        double fitness_squared;
        for(int i = 0; i < popsize; i++)
        {
            if (population[i] > 0)
            {
                fitness = fitness_calculator(population[i]);
            }
            else
            {
                fitness = fitness_calculator(population[i] * -1) ;
            }
            fitness_squared = pow(fitness,2);
            sum += fitness;
            sum_squared += fitness_squared;
        }
        return ((sum_squared) - (pow(sum,2)) / popsize) /popsize;
    }

    // calculates fitness variance of population
    double fitness_category_variance_calculator()
    {
        double sum = 0;
        double sum_squared = 0;
        for(int i = 0; i < popsize; i++)
        {
            sum += population[i];
            sum_squared += pow(population[i],2);
        }
        return ((sum_squared) - (pow(sum,2)) / popsize) /popsize;
    }


    /** Function that runs a neutral theory on "population" and saves results at intervals**/

    // runs simulation for seed, generation, strenght, mutation, rate population, pop interval, can specify directory path in which to ave.
    void runsimulation(long generations, long mean_interval, long variance_interval, long abundances_interval)
    {
        // build recycler
        reset_recycler();
        if(type_abundances.size() != type_categories.size()){cout << "pre sim error 1" << endl;}
        if(types.size() != population.size()){cout << "pre sim error 2" << endl;}
        if(population.size() != (unsigned)std::accumulate(type_abundances.begin(), type_abundances.end(), 0)){cout << "pre sim error 3" << endl;}
        for(unsigned int i = 0; i < types.size(); i++)if(types[i] < 0){cout <<"pre sim error 4" << endl;}

        popsize = population.size();
        min_cat = *min_element(type_categories.begin(), type_categories.end());
        max_cat = *max_element(type_categories.begin(), type_categories.end());

        fitness_category_means.push_back(fitness_category_mean_calculator());
        fitness_means.push_back(fitness_mean_calculator());
        fitness_variance.push_back(fitness_variance_calculator());
        fitness_category_variance.push_back(fitness_category_variance_calculator());
        vector_of_type_abundances.push_back(type_abundances);
        vector_of_type_categories.push_back(type_categories);
        extinct_counter.push_back(extinct_type_counter);


        // runs simulation update outputs.
        for( int a = 1; a <= generations; a = a + 1 )
        {
            generationstep();
             // reset recyler (to avoid memory overload);
            reset_recycler();



            if (a % mean_interval == 0)
            {

                fitness_means.push_back(fitness_mean_calculator());
                fitness_category_means.push_back(fitness_category_mean_calculator());
            }
            if (a % variance_interval == 0)
            {
                fitness_variance.push_back(fitness_variance_calculator());
                fitness_category_variance.push_back(fitness_category_variance_calculator());
            }
            if (a % abundances_interval == 0)
            {
                vector_of_type_abundances.push_back(type_abundances);
                vector_of_type_categories.push_back(type_categories);
                extinct_counter.push_back(extinct_type_counter);

            }
        }
        if(type_abundances.size() != type_categories.size()){cout << "post sim error 1" << endl;}
        if(types.size() != population.size()){cout << "post sim error 2" << endl;}
        if(population.size() != (unsigned)std::accumulate(type_abundances.begin(), type_abundances.end(), 0)){cout << "post sim error 3" << endl;}
        for(unsigned int i = 0; i < types.size(); i++)if(types[i] < 0){cout <<"post sim error 4" << endl;}

        // write
        string fitness_mean_path = temporary_DirectoryPath + "Fitness_mean_" + to_string(seed) + "_" + to_string(task_number) +".csv";; // define path
        string fitness_variance_path = temporary_DirectoryPath + "Fitness_Var_" + to_string(seed) + "_" + to_string(task_number) +".csv"; // define path
        string fitness_category_mean_path = temporary_DirectoryPath + "Cat_mean_" + to_string(seed) + "_" + to_string(task_number) +".csv"; // define path
        string fitness_category_variance_path = temporary_DirectoryPath + "Cat_Var_" + to_string(seed) + "_" + to_string(task_number) +".csv"; // define path
        string abundances_path = temporary_DirectoryPath + "Type_abund_" + to_string(seed) + "_" + to_string(task_number) +".csv"; // define path
        string categories_path = temporary_DirectoryPath + "Type_Cat_" + to_string(seed) + "_" + to_string(task_number) +".csv";
        string extinct_path = temporary_DirectoryPath + "Type_extinct_" + to_string(seed) + "_" + to_string(task_number) +".csv";


        // and save everything using dynamic system
        vectors_of_longs =  vector_of_type_abundances;
        savevectorlongs(abundances_path);
        vectors_of_longs.clear();


        vectors_of_longs =  vector_of_type_categories;
        savevectorlongs(categories_path);
        vectors_of_longs.clear();

        longs = extinct_counter;
        savelongs(extinct_path);
        longs.clear();

        doubles = fitness_means;
        savedoubles(fitness_mean_path);
        doubles.clear();

        doubles = fitness_variance;
        savedoubles(fitness_variance_path);
        doubles.clear();


        doubles = fitness_category_variance;
        savedoubles(fitness_category_variance_path);
        doubles.clear();

        doubles = fitness_category_means;
        savedoubles(fitness_category_mean_path);
        doubles.clear();

    }
    /** Function that sets initail values for private objects**/

    // setter, sets all the private variables.
    void setup(int internal_seed_number, int internal_task_number,double internal_selection_strength, double internal_mutation_rate, long internal_home_popsize)
    {

        //set paramaters that are parsed
        seed = internal_seed_number;
        task_number = internal_task_number;
        selection_strength =internal_selection_strength;
        mutation_rate = internal_mutation_rate;

        // home
        home_popsize = internal_home_popsize ;
        home_generations = home_popsize*4; // run for popsize number of generations for burnin
        home_interval_mean = home_popsize/100;
        home_interval_variance = home_popsize/100;
        home_interval_abundances = home_popsize/10;


        // assign default values to populations
        home_population.assign(home_popsize,100);
        home_types.assign(home_popsize,0);
        home_type_abundances.assign(1,home_popsize);
        home_type_categories.assign(1,100);

        // set up temporary objects
        population = home_population;
        types = home_types;
        type_abundances = home_type_abundances;
        type_categories = home_type_categories;
                // paths
       DirectoryPath = "DFE1_Basic_" ;

        check_max_cat =false ;// default
        check_min_cat = false;// default
        extinct_type_counter = 0;
        // set seed
        rng.set_seed(seed);


    }

    /** Function thats called to piece verything together**/

    void run(int seed,int task_number, double internal_selection_strength, double internal_mutation_rate,long internal_home_popsize)
    {
        // set up times

        setup(seed, task_number,internal_selection_strength, internal_mutation_rate, internal_home_popsize);

        /** burn in**/
        temporary_DirectoryPath = DirectoryPath;

//

        // run burn in
        cout << "Burnin Begin" << endl;
        runsimulation(home_generations,home_interval_mean,home_interval_variance,home_interval_abundances);
        saveparams();
        cout << "Burnin End" << endl;

        /** invador evolution**/
        //clear
        home_population.clear();
        home_types.clear();
        home_type_abundances.clear();
        home_type_categories.clear();

        // clear everything at the end
        population.clear();
        types.clear();
        type_abundances.clear();
        type_categories.clear();

        // to save
        fitness_means.clear();
        fitness_category_means.clear();
        fitness_variance.clear();
        fitness_category_variance.clear();
        vector_of_type_abundances.clear();
        vector_of_type_categories.clear();
        extinct_counter.clear();

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
    long seed_number = -1; // the task number
    double wall_time = -1; // the max time (minutes)

//    // read command line data
    seed_number = jobconvertor(argv[1]);
    wall_time = jobconvertor(argv[2])  * 3600; //(convert hours to seconds) ;

    vector<long> jm_vect;
    vector<double> nu_vect; // mutation rates
    vector<double> s_vect; // selection strength

    jm_vect.clear();
    nu_vect.clear();
    s_vect.clear();

    // uninteresting paramater

    jm_vect.push_back(500);
    jm_vect.push_back(1000);
    jm_vect.push_back(2500);
    jm_vect.push_back(5000);
    jm_vect.push_back(7500);
    jm_vect.push_back(10000);
    jm_vect.push_back(25000);
    jm_vect.push_back(50000);
    jm_vect.push_back(75000);
    jm_vect.push_back(100000);

    nu_vect.push_back(0);
    nu_vect.push_back(0.00001);
    nu_vect.push_back(0.00005);
    nu_vect.push_back(0.0001);
    nu_vect.push_back(0.00025);
    nu_vect.push_back(0.0005);
    nu_vect.push_back(0.00075);
    nu_vect.push_back(0.001);
    nu_vect.push_back(0.0025);
    nu_vect.push_back(0.005);
    nu_vect.push_back(0.0075);
    nu_vect.push_back(0.01);
    nu_vect.push_back(0.025);
    nu_vect.push_back(0.05);
    nu_vect.push_back(0.075);
    nu_vect.push_back(0.1);
    nu_vect.push_back(0.25);
    nu_vect.push_back(0.5);
    nu_vect.push_back(0.75);
    nu_vect.push_back(1);

    s_vect.push_back(0.1);
    s_vect.push_back(0.05);
    s_vect.push_back(0.01);
    s_vect.push_back(0.005);
    s_vect.push_back(0.001);
    s_vect.push_back(0.0005);
    s_vect.push_back(0.0001);
    s_vect.push_back(0.00005);
    s_vect.push_back(0.00001);
    s_vect.push_back(0);

    vector<long> jm_list;
    vector<double> nu_list;
    vector<double> s_list;

    // construct paramater combinations for selection and mutation

    for (unsigned int k = 0 ; k < s_vect.size() ; k ++)
        {
            for (unsigned int j = 0 ; j < nu_vect.size() ; j ++)
            {

//            for (unsigned int i = 0 ; i < jm_vect.size() ; i ++)
//            {
//                jm_list.push_back(jm_vect[i]);
                nu_list.push_back(nu_vect[j]);
                s_list.push_back(s_vect[k]);
//            }
        }
    }

    cout << seed_number << " = seed \n";
    cout << wall_time << " = max time \n";

    Simulation simulation;


    // keep running until time runs out
    time_t start , end;
    time(&start);
    bool timeout = false;
    int task_number = 1;
   do{
        for (unsigned int i = 0 ; i < jm_vect.size() ; i ++)
        {
            cout << "starting task " << task_number << " with popsize " << jm_vect[i] << endl;
            simulation.run(seed_number,task_number,s_list[seed_number-1],nu_list[seed_number-1],jm_vect[i]);
            time(&end);
            cout << "Task took " <<  difftime(end,start) << "second" << endl;
            task_number++;
        if (difftime(end,start) >= (wall_time)){timeout = true;}
        if(timeout){break;}
        }
    if(task_number >= 200){break;}
    }
    while (!timeout);

    cout << "total time is " <<  difftime(end,start) << "second" << endl;
    return 0;
}

