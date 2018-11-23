/* Indet_growth_2exp.h -    Header file specifying the elementary life-history functions of
                            a structured consumer-unstructured resource model, in which the
                            consumer follows a net-production DEB model with arbitrary, allometric
                            scaling functions of ingestion and maintenance as a function of body size.
 
                            Functional Response is Holling Type II
 
                            Layout for the i-state variables:
 
                                istate[i][0] : Age(i)
                                istate[i][1] : Size(i)
 
                            Layout for the environmental variables:
 
                                E[0] : Resource
 
                            Layout for the interaction (and output) variables:
 
                                I[0][0] : Total ingestion
                                I[0][1] : Juvenile biomass
                                I[0][2] : Adult biomass (growing and non-growing)
                                I[0][3] : Total Biomass
                                I[0][4] : Age at Maturation * Birth Rate
                                I[0][5] : Asymptotic Size * Birth Rate - 1
                                I[0][6] : Population Level Reproduction Rate
 
                            Last modification: VH - October 4, 2015
 */

// Which allocation function to use. 0 for linear, 1 for sigmoid and 2 for exponential
#define ALLO_FUNC   1

#if (ALLO_FUNC == 0)
#warning ALLOCATION FUNCTION IS LINEAR
#elif (ALLO_FUNC == 1)
#warning ALLOCATION FUNCTION IS SIGMOID
#elif (ALLO_FUNC == 2)
#warning ALLOCATION FUNCTION IS EXPONENTIAL
#endif

#define MORT_SIZE	0

#if (MORT_SIZE == 0)
#warning MORTALITY IS CONSTANT
#elif (MORT_SIZE == 1)
#warning MORTALITY IS SIZE-DEPENDENT
#endif

/*
 *===========================================================================
 * 		DEFINITION OF PROBLEM DIMENSIONS AND NUMERICAL SETTINGS
 *===========================================================================
 */
// Dimension settings: Required
#define POPULATION_NR		1
#define STAGES              2
#define	I_STATE_DIM         2
#define	ENVIRON_DIM         1
#define INTERACT_DIM		9
#define	PARAMETER_NR		15

// Numerical settings: Optional (default values adopted otherwise)
#define MIN_SURVIVAL		1.0E-9		// Survival at which individual is considered dead
#define MAX_AGE             100000		// Give some absolute maximum for individual age

#define DYTOL               1.0E-7		// Variable tolerance
#define RHSTOL              1.0E-6		// Function tolerance
#define ESSTOL				1.0E-6		// ESS tolerance

#define COHORT_NR           1000
#define ALLOWNEGATIVE       0

/*
 *===========================================================================
 * 		DEFINITION OF ALIASES
 *===========================================================================
 */
// Define aliases for the istate variables
#define AGE         istate[0][0]
#define SIZE		istate[0][1]

// Define aliases for the environmental variables
#define R           E[0]

// Define aliases for the parameters
#define RMAX		parameter[ 0]	// Default: 30      Maximum Resource Density
#define DELTA		parameter[ 1]	// Default: 0.1     Resource turn-over rate

#define M           parameter[ 2]	// Default: 0.1		Maximum Ingestion rate of adult (s = sj)
#define Q           parameter[ 3]	// Default: 1.0		Maximum ingestion exponent

#define T           parameter[ 4]	// Default: 0.01	Maintenance rate of adult (s = sj)
#define P           parameter[ 5]	// Default: 1		Maintenance rate exponent

#define MU   		parameter[ 6]   // Default: 0.0015	Mortality rate of adult (s = sj)
#define MUJ   		parameter[ 7]   // Default: 0.0     Mortality rate juveniles
#define MUA   		parameter[ 8]   // Default: 0.0     Mortality rate adults

#define XB          parameter[ 9]	// Default: 0.1 	Mass of newborn individual (gram)
#define XJ          parameter[10]   // Default: 1       Mass of adult
#define Xref        parameter[11]   // Default: 1       Scaling reference size
#define XM          parameter[12]	// Default: 10		Maximum size of adult (gram)

#define SIGMA		parameter[13]	// Default: 0.5		Food conversion efficiency
#define H           parameter[14]   // Default: 3.0     Half-saturation constant

/*
 *===========================================================================
 * 		DEFINITION OF NAMES AND DEFAULT VALUES OF THE PARAMETERS
 *===========================================================================
 */
// At least two parameters should be specified in this array
char  *parameternames[PARAMETER_NR] =
{"RMAX", "DELTA", "M", "Q", "T", "P", "MU", "MUJ", "MUA", "XB", "XJ", "Xref", "XM", "SIGMA", "H"};

// These are the default parameters values
double	parameter[PARAMETER_NR] =
{30, 0.1, 1.0, 1.0, 0.1, 1.0, 0.015, 0, 0, 0.1, 1.0, 1.0, 10.0, 0.5, 3.0};

/*
 *===========================================================================
 * 		DEFINITION OF THE LIFE HISTORY MODELS FOLLOWS BELOW
 *===========================================================================
 * Specify the number of states at birth for the individuals in all structured
 * populations in the problem in the vector BirthStates[].
 *===========================================================================
 */

void SetBirthStates(int BirthStates[POPULATION_NR], double E[])
{
    BirthStates[0] = 1;
    
    return;
}

/*
 *===========================================================================
 * Specify all the possible states at birth for all individuals in all
 * structured populations in the problem. BirthStateNr represents the index of
 * the state of birth to be specified. Each state at birth should be a single,
 * constant value for each i-state variable.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 *===========================================================================
 */

void StateAtBirth(double *istate[POPULATION_NR], int BirthStateNr, double E[])
{
    
    AGE = 0.0;
    SIZE = XB;
    
    return;
}


/*
 *===========================================================================
 * Specify the threshold determining the end point of each discrete life
 * stage in individual life history as function of the i-state variables and
 * the individual's state at birth for all populations in every life stage.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 *===========================================================================
 */

void IntervalLimit(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
                   double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
                   double limit[POPULATION_NR])
{
    switch (lifestage[0])
    {
        case 0:
            limit[0] = SIZE - XJ;
            break;
    }
    
    return;
}




/*
 *===========================================================================
 * Specify the individual development of individuals as function of i-state
 * and environment for all individuals of all populations in every life stage
 *
 * Notice that the first index of the variables 'istate[][]' and 'growth[][]'
 * refers to the number of the structured population, the second index refers
 * to the number of the individual state variable. The interpretation of the
 * latter is up to the user.
 *===========================================================================
 */

void Development(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
            double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
            double growth[POPULATION_NR][I_STATE_DIM])
{
    double	ingest, netproduction, netproduction_plus;
    double	kappa;
    
#if (ALLO_FUNC == 0)
    kappa       = (XM - SIZE)/(XM - XJ);
#elif (ALLO_FUNC == 1)
    double      L;
    L           = (SIZE - XJ)/(XM - XJ);
    kappa       = 1 - 3*L*L + 2*L*L*L;
#elif (ALLO_FUNC == 2)
    double      L2, L1;
    L2          = log(10000) / (XM - XJ);
    L1          = 1 / exp(-L2);
    kappa       = L1 * exp(-L2*SIZE);
#endif
    
    growth[0][0]        = 1.0;
    
    ingest              = M*pow(SIZE / Xref, Q) * R/(H+R);
    netproduction       = SIGMA*ingest - T*pow(SIZE / Xref, P);
    netproduction_plus  = max(netproduction, 0.0);
    
    if (lifestage[0] == 0)
    {
        growth[0][1]    = netproduction_plus;
    }
    else
    {
        growth[0][1]    = kappa*netproduction_plus;
    }
    
    return;
}


/*
 *===========================================================================
 * Specify the possible discrete changes (jumps) in the individual state
 * variables when ENTERING the stage specified by 'lifestage[]'.
 *
 * Notice that the first index of the variables 'istate[][]' and 'growth[][]'
 * refers to the number of the structured population, the second index refers
 * to the number of the individual state variable. The interpretation of the
 * latter is up to the user.
 *===========================================================================
 */

void DiscreteChanges(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
                     double *birthstate[POPULATION_NR], int BirthStateNr, double E[])
{
    return;
}


/*
 *===========================================================================
 * Specify the fecundity of individuals as a function of the i-state
 * variables and the individual's state at birth for all populations in every
 * life stage.
 *
 * The number of offspring produced has to be specified for every possible
 * state at birth in the variable 'fecundity[][]'. The first index of this
 * variable refers to the number of the structured population, the second
 * index refers to the number of the birth state.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 *===========================================================================
 */

void Fecundity(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
               double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
               double *fecundity[POPULATION_NR])
{
    double	ingest, netproduction, netproduction_plus;
    double	kappa;
    
#if (ALLO_FUNC == 0)
    kappa       = (XM - SIZE)/(XM - XJ);
#elif (ALLO_FUNC == 1)
    double      L;
    L           = (SIZE - XJ)/(XM - XJ);
    kappa       = 1 - 3*L*L + 2*L*L*L;
#elif (ALLO_FUNC == 2)
    double      L2, L1;
    L2          = log(10000) / (XM - XJ);
    L1          = 1 / exp(-L2);
    kappa       = L1 * exp(-L2*SIZE);
#endif
    
    ingest              = M*pow(SIZE / Xref, Q) * R/(H+R);
    netproduction       = SIGMA*ingest - T*pow(SIZE / Xref, P);
    netproduction_plus  = max(netproduction, 0.0);

    
    if (lifestage[0] == 0)
    {
        fecundity[0][0]    = 0.0;
    }
    else
    {
        fecundity[0][0]    = (1-kappa)*netproduction_plus / XB;
    }
    
    return;
}



/*
 *===========================================================================
 * Specify the mortality of individuals as a function of the i-state
 * variables and the individual's state at birth for all populations in every
 * life stage.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 *===========================================================================
 */

void Mortality(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
               double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
               double mortality[POPULATION_NR])
{
	
// 	double	ingest, netproduction, netproduction_min;
//     double	kappa;
//     
// #if (ALLO_FUNC == 0)
//     kappa       = (XM - SIZE)/(XM - XJ);
// #elif (ALLO_FUNC == 1)
//     double      L;
//     L           = (SIZE - XJ)/(XM - XJ);
//     kappa       = 1 - 3*L*L + 2*L*L*L;
// #elif (ALLO_FUNC == 2)
//     double      L2, L1;
//     L2          = log(10000) / (XM - XJ);
//     L1          = 1 / exp(-L2);
//     kappa       = L1 * exp(-L2*SIZE);
// #endif
//     
//     ingest              = M*pow(SIZE / XJ, Q) * R/(H+R);
//     netproduction       = SIGMA*ingest - T*pow(SIZE / XJ, P);
//     netproduction_min	= min(netproduction, 0.0);
	
#if (MORT_SIZE == 0)
    if (lifestage[0] == 0)
    {
       	mortality[0] = MU + MUJ;
    }
    else
    {
       	mortality[0] = MU + MUA;
    }
#elif (MORT_SIZE == 1)	
		mortality[0] = MU + MUJ*exp(-SIZE / XJ) + MUA / (1 + XM * exp(-SIZE));
#endif
	
    return;
}


/*
 *===========================================================================
 * For all the integrals (measures) that occur in interactions of the
 * structured populations with their environments and for all the integrals
 * that should be computed for output purposes (e.g. total juvenile or adult
 * biomass), specify appropriate weighing function dependent on the i-state
 * variables, the environment variables and the current life stage of the
 * individuals. These weighing functions should be specified for all
 * structured populations in the problem. The number of weighing functions
 * is the same for all of them.
 *
 * Notice that the first index of the variables 'istate[][]' and 'impact[][]'
 * refers to the number of the structured population, the second index of the
 * variable 'istate[][]' refers to the number of the individual state variable,
 * while the second index of the variable 'impact[][]' refers to the number of
 * the interaction variable. The interpretation of these second indices is up
 * to the user.
 *===========================================================================
 */

void Impact(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
            double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
            double impact[POPULATION_NR][INTERACT_DIM])
{
    double	ingest, netproduction, netproduction_plus;
    double	kappa;
    double  mort_j;
    
#if (ALLO_FUNC == 0)
    kappa       = (XM - SIZE)/(XM - XJ);
#elif (ALLO_FUNC == 1)
    double      L;
    L           = (SIZE - XJ)/(XM - XJ);
    kappa       = 1 - 3*L*L + 2*L*L*L;
#elif (ALLO_FUNC == 2)
    double      L2, L1;
    L2          = log(10000) / (XM - XJ);
    L1          = 1 / exp(-L2);
    kappa       = L1 * exp(-L2*SIZE);               
#endif

#if (MORT_SIZE == 0)
       	mort_j = MU + MUJ;
#elif (MORT_SIZE == 1)	
		mort_j = MU + MUJ*exp(-SIZE / XJ) + MUA / (1 + XM * exp(-SIZE));
#endif
    
    ingest              = M*pow(SIZE / Xref, Q) * R / (H + R);
    netproduction       = SIGMA*ingest - T*pow(SIZE / Xref, P);
    netproduction_plus  = max(netproduction, 0.0);

    switch (lifestage[0])
    {
        case 0:
            impact[0][0] = ingest;                               // Ingestion
            impact[0][1] = SIZE;                                 // Juvenile Biomass
            impact[0][2] = 0.0;                                  // Adult Biomass
            impact[0][3] = SIZE;                                 // Total Biomass
            impact[0][4] = 0.0;                                  // Asymptotic Size
            impact[0][5] = 0.0;                                  // Reproduction
            impact[0][6] = netproduction_plus - mort_j*SIZE;     // Total Juvenile Production - Total Juvenile Mortality
            impact[0][7] = mort_j / Survival(0);                  // ln(- Juvenile Survival * b)
            impact[0][8] = 0.0;                       // Age at Maturation * b            
           
            break;
        case 1:
            impact[0][0] = ingest;                                  // Ingestion
            impact[0][1] = 0.0;                                     // Juvenile Biomass
            impact[0][2] = SIZE;                                 	// Adult Biomass
            impact[0][3] = SIZE;                                 	// Total Biomass
            impact[0][4] = (kappa * netproduction_plus / Survival(0));  // Integrates Growth Function. Gives Asymptotic Size * Birth Rate - XJ
            impact[0][5] = ((1-kappa) * netproduction_plus);        // Population Level Reproduction in biomass
            impact[0][6] = 0.0;                                     // Total Juvenile Production - Total Juvenile Mortality
            impact[0][7] = 0.0;                                     // Juvenile Survival
            impact[0][8] = ((1-kappa) * netproduction_plus / XB);	// b * Individual level fecundity in offspring numbers

            break;
          /* In equilibrium the following equilities holds; 
             *
             *  Population Level Maturation - Population Level Reproduction = Total Juvenile Net Production - Total Juvenile Mortality      
             *  I[0][6] calculates the RHS.
             *  I[0][5] calculates the population level reproduction
             *  Population Level Maturation follows from I[0][6] + I[0][5]
             *
             *  I[0][7] calculates Juvenile Survival * birthrate
             *  I[0][8] calculates Age at Maturation * birthrate
             */
    }
    return;
}


/*
 *===========================================================================
 * Specify the type of each of the environment variables by setting
 * the entries in EnvironmentType[ENVIRON_DIM] to PERCAPITARATE, GENERALODE
 * or POPULATIONINTEGRAL based on the classification below:
 *
 * Set an entry to PERCAPITARATE if the dynamics of E[j] follow an ODE and 0
 * is a possible equilibrium state of E[j]. The ODE is then of the form
 * dE[j]/dt = P(E,I)*E[j], with P(E,I) the per capita growth rate of E[j].
 * Specify the equilibrium condition as condition[j] = P(E,I), do not include
 * the multiplication with E[j] to allow for detecting and continuing the
 * transcritical bifurcation between the trivial and non-trivial equilibrium.
 *
 * Set an entry to GENERALODE if the dynamics of E[j] follow an ODE and 0 is
 * NOT an equilibrium state of E. The ODE then has a form dE[j]/dt = G(E,I).
 * Specify the equilibrium condition as condition[j] = G(E,I).
 *
 * Set an entry to POPULATIONINTEGRAL if E[j] is a (weighted) integral of the
 * population distribution, representing for example the total population
 * biomass. E[j] then can be expressed as E[j] = I[p][i]. Specify the
 * equilibrium condition in this case as condition[j] = I[p][i].
 *
 * Notice that the first index of the variable 'I[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the interaction variable. The interpretation of the latter
 * is up to the user. Also notice that the variable 'condition[j]' should
 * specify the equilibrium condition of environment variable 'E[j]'.
 *===========================================================================
 */

const int EnvironmentType[ENVIRON_DIM] = {GENERALODE};

void EnvEqui(double E[], double I[POPULATION_NR][INTERACT_DIM],
             double condition[ENVIRON_DIM])
{
    condition[0] = DELTA*(RMAX - R) - I[0][0];
    
    return;
}

/*==============================================================================*/

