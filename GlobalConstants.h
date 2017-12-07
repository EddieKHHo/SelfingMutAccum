//Global constants for the simulation


extern std::string   MODELFILENAME;

extern double   NUMOFREPS;
extern double   NUMOFGENERATIONS;
extern double   MAXNUMOFGENERATIONS;

extern int   	L;
extern int		POPULATION;
extern double   MUTATIONRATE;
extern double   H;
extern double   S;
extern double   K;

extern double   INITIALSELFING;
extern double   SELFINGMUTATIONRATE;
extern double   SIGMAS;
extern double   FINALSELFINGMUTATIONRATE;
extern double   FINALSIGMAS;

extern double   INITIALRECOMBINATION;
extern double   RECOMBMUTATIONRATE;
extern double   SIGMAR;

//**********Eddie MODDED Dec 30 2014; Apr 24 2015
extern double   STOCHASTIC_FLUCTUATION;

//Deterministic fluctuations
extern double   PERIOD;
extern double   T0;
extern double   T1;
extern double   TAU;

//Stochastic fluctuations
extern double	F_AC;
extern double	PHI;
extern double	PSI;
extern double	EPSILON;


//variable to indicate the frequency we want to save main data
extern int      SAVE_DATA_INTERVAL;


//**********Eddie MODDED July 21 2015
extern int		RUNTIME_ERRORS;
extern double 	P_LOCI_GROUP[3];

extern double	DELTA_S0_A;
extern double	DELTA_S0_B;
extern double	DELTA_S1_A;
extern double	DELTA_S1_B;
extern double	DELTA_S2_C;

//Proportional loci constant selected
extern double	P_CONSTANT;
//Proportional of fluctuating selected loci in mut group 0 (rather than group 1)
extern double	P_FLUC_G0;
extern double	ENV_FLUC_INDICATOR;
//**********Eddie MODDED July 21 2015

//**********Eddie MODDED Mar 8 2016
extern double CHECK_FIXATION;
//**********Eddie MODDED Mar 8 2016

//**********Eddie MODDED Mar 8 2016
extern double SAVE_DATA_MIN_GEN;
extern double Q_NEUTRAL_LOCI;
//**********Eddie MODDED Mar 8 2016

//**********Eddie MODDED May 3 2016
extern double ETA;
//**********Eddie MODDED May 3 2016

extern double SELF_LENGTH;
extern double IGNORE_LOCI;

extern double PFIX_MODEL;
extern double PFIX_H;
extern double PFIX_S;
extern double PFIX_K;

extern double U_BEN;
extern double S_BEN;
extern double H_BEN;

#define PROB_BEN 0.0