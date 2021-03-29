# pertussis_main.py  
# Andrés Escallón
# Marzo 2021

#from escbox.py import *

from pertussis_header import *
from pertussis_test import *

'''
/*
 *==========================================================================
 *
 *		LABELLING ENVIRONMENT AND I-STATE VARIABLES
 *
 *==========================================================================
 */
 '''

time = env[0]
S =	env[1]
R = env[2]
age = i_state(0)


'''
/*
 *==========================================================================
 *
 *		DEFINING AND LABELLING CONSTANTS AND PARAMETERS
 *
 *==========================================================================
 */
 '''

MU = parameter[0]
N = parameter[1]
R0 = parameter[2]
INF_PERIOD = parameter[3]
LAT_PERIOD = parameter[4]


alpha=0


'''
/*
 *==========================================================================
 *
 * USER INITIALIZATION ROUTINE ALLOWS OPERATIONS ON INITIAL POPULATIONS
 *
 *==========================================================================
 */
'''
#def	UserInit( int argc, char **argv, double *env,  population *pop)
def UserInit(argc, argv, env, pop):

  I = 0.0
  for i in range(cohort_no[0]): I += pop[0][i][number]
  S = N - I - R

  alpha = R0*MU*MU*MU*INF_PERIOD/N/(2*exp(-MU*(LAT_PERIOD+INF_PERIOD))+exp(-MU*(LAT_PERIOD+INF_PERIOD))*MU*INF_PERIOD-2*exp(-MU*LAT_PERIOD)+exp(-MU*LAT_PERIOD)*MU*INF_PERIOD)

  



'''
/*
 *==========================================================================
 *
 *	SPECIFICATION OF THE NUMBER AND VALUES OF BOUNDARY POINTS
 *
 *==========================================================================
 */
'''
#void	SetBpointNo(double *env, population *pop, int *bpoint_no)
def SetBpointNo(env, pop, bpoint_no):
  
  bpoint_no[0]=1
  

#/*==========================================================================*/

#void	SetBpoints(double *env, population *pop, population *bpoints)
def SetBpoints(env, pop, bpoints):

  bpoints[0][0][age]=0.0
  
'''
/*
 *==========================================================================
 *
 *			SPECIFICATION OF DERIVATIVES
 *
 *==========================================================================
 */
'''

#void	Gradient(double *env,     population *pop,     population *ofs,
#		 double *envgrad, population *popgrad, population *ofsgrad,
#		 population *bpoints)

def Gradient(env, pop, ofs, envgrad, popgrad, ofsgrad, bpoints):


  I = 0.0
  InfectionForce = 0.0
  
  for i in range(cohort_no[0]):
    
      I  += pop[0][i][number]
      tau = pop[0][i][age] - LAT_PERIOD

      if (tau > 0.0):
	    tmp = alpha*tau*(1.0 - tau/INF_PERIOD)*pop[0][i][number]
      else:
	    tmp = 0.0

      if (tmp > 0.0): InfectionForce += tmp

      popgrad[0][i][number] = -MU*pop[0][i][number]
      popgrad[0][i][age]    = 1
    

  if(ofs[0][0][number] > 1.0E-10):#		/* cohort if non-zero       */
    
      I += ofs[0][0][number]
      tau = ofs[0][0][age]/ofs[0][0][number] - LAT_PERIOD

      if (tau > 0.0):
	    tmp = alpha*tau*(1.0 - tau/INF_PERIOD)*ofs[0][0][number]
      else:
	    tmp = 0.0

      if (tmp > 0.0): InfectionForce += tmp
    

  envgrad[0] = 1
  envgrad[1] = MU*(I + R) - InfectionForce*S
  envgrad[2] = -MU*R
  
						#/* The derivatives for the  */
						#/* boundary cohort	    */
  ofsgrad[0][0][number] = -MU*ofs[0][0][number] + InfectionForce*S
  ofsgrad[0][0][age]    = -MU*ofs[0][0][age] + ofs[0][0][number]

  



'''
/*
 *==========================================================================
 *
 *	SPECIFICATION OF EVENT LOCATION AND DYNAMIC COHORT CLOSURE
 *
 *==========================================================================
 */
'''

#void	EventLocation(double *env, population *pop, population *ofs,
#		      population *bpoints, double *events)
def EventLocation(env, pop, ofs, bpoints, events):
    print ("EventLocation")
  

#/*==============================================================================*/

#int	ForceCohortEnd(double *env, population *pop, population *ofs,
#		       population *bpoints)
def ForceCohortEnd(env, pop, ofs, bpoints):

  return NO_COHORT_END

'''
/*
 *==========================================================================
 *
 *		SPECIFICATION OF BETWEEN COHORT CYCLE DYNAMICS
 *
 *==========================================================================
 */
 '''

#void	InstantDynamics(double *env, population *pop, population *ofs)
def	InstantDynamics(env, pop, ofs):

   
  for i in range(cohort_no[0]):#			/* Empty cohorts that have  */
    						#/* reached maximum age      */
      if(pop[0][i][age] > (LAT_PERIOD + INF_PERIOD)):
	
	    R += pop[0][i][number]
	    pop[0][i][number] = -1.0
	
'''
/*
 *==========================================================================
 *
 *			SPECIFICATION OF OUTPUT VARIABLES
 *
 *==========================================================================
 */
'''

#void	DefineOutput(double *env, population *pop, double *output)
def DefineOutput(env, pop, output):
  
  E = 0.0
  I = 0.0
  InfectionForce = 0.0
  
  for i in range(cohort_no[0]):
    
      tau = pop[0][i][age] - LAT_PERIOD

      if (tau > 0.0):
	
	    tmp = alpha*tau*(1.0 - tau/INF_PERIOD)*pop[0][i][number]
	    I += pop[0][i][number]
	  if (tmp > 0.0): InfectionForce += tmp
	
      else: E += pop[0][i][number]
    

  output[0] = env[1]
  output[1] = I
  output[2] = R
  output[3] = E
  output[4] = InfectionForce
  output[5] = S + I + E + R

  




#/*==========================================================================*/
