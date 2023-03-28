/*********************************************
 * OPL 22.1.0.0 Model
 * Author: emazied
 * Creation Date: Dec 5, 2022 at 7:46:27 PM
 *********************************************/

//// Predefined Scalar Data parameters
int P = ...;//The cardinality of CPU resource set ... max number of available cpu cores, i.e., 64 CPU cores

int G = ...;//The cardinality of GPU resource set ... max number of available gpu cores, i.e., 8  GPU cores 

int a_p = ...;//The memory page size of the memory that is available for CPU architecture, ... Bytes

int a_g = ...;//The memory page size of the memory that is available for GPU architecture, ... Bytes

int M_p = ...;//Max number of available memory pages for CPU ops, i.e., ... pages

int M_g = ...;//Max number of available memory pages for GPU ops, i.e., ... pages

int I = ...;//Cardinality of containers set, i.e., 3 containers {enumerate > 1: eMBB service, 2: mMTC service, 3: urllc service}

int z = ...;//Data size,i.e.,since we assume matrix operations on binary data elements (1 or 0), we consider z = 4 Bytes

float pi_p = ...;   //Peak performance of CPU architecture, ... Operations/second

float pi_g = ...;   //Peak performance of GPU architecture, ... Operations/second

float beta_p = ...; //Memory bandwidth of CPU architecture, ... Bytes/second

float beta_g = ...; //Memory bandwidth of GPU architecture, ... Bytes/second

float erf_inv = ...; // error inverse function of alpha


//// Ranges ////

range cpu_cores = 1..P;

range gpu_cores = 1..G;


range cpu_mem_pages = 1..M_p;

range gpu_mem_pages = 1..M_g;

range service_categories = 1..I;

////////////////////////////////////////////////////////////////////////////////////////

//// vectors //////

// cost vector that represent the cost of using MEC resource (p, g, m_p, or m_g) during a period (i.e., scaling period) unit cost for a container that supports service category

float cpu_cost[service_categories][cpu_cores];

float gpu_cost[service_categories][gpu_cores];

float cpu_mem_cost[service_categories][cpu_mem_pages];

float gpu_mem_cost[service_categories][gpu_mem_pages];

////////////////////////////////////////////////////////////////////////

/////////// vector of the parameters that reflect stochastic behavior of the model

int h[service_categories] = ...; 
int n[service_categories] = ...;

float gamma[service_categories] = ...;
float rho[service_categories] = ...;

float ops[service_categories]; // number of running operations per container

float delta[service_categories] = ...;

float avg_r[service_categories] = ...;

float var_r[service_categories] = ...;

float phi_r[service_categories];


///////////////////////////////////////////////////////////

///// Parameters for processing time constraint operations 
float u[service_categories];
float v_p[service_categories];
float v_g[service_categories];
float M1[service_categories];

//// to compute the resource limits ....
int cpu_limit[service_categories];
int memc_limit[service_categories];

int gpu_limit[service_categories];
int memg_limit[service_categories];

/// Routine for computing the ubber and lower bounds of process time constraints process_time_constraint_lb[service_categories], and process_time_constraint_ub[service_categories]
execute
{  
for (var i in service_categories) 
{
  ops[i] = 2*h[i]*rho[i]+4*n[i]*gamma[i]+n[i]; 
  phi_r[i] = ((Opl.sqrt(2*var_r[i]))*erf_inv) + (avg_r[i]);


  u[i] =  avg_r[i] + phi_r[i]*Opl.sqrt(var_r[i]);
    
  v_p[i] = Opl.minl(beta_p*u[i]*ops[i]/(h[i]*n[i]*z), pi_p);
  v_g[i] = Opl.minl(beta_g*u[i]*ops[i]/(h[i]*n[i]*z), pi_g);
  

M1[i] = Opl.maxl(Opl.maxl( (Opl.maxl(G, u[i]*ops[i]-delta[i]*v_g[i]*G) ), (Opl.maxl(P, u[i]*ops[i]-delta[i]*v_p[i]*P) )), u[i]*ops[i]);

  
} // end of for var i scope
 
for (var g in gpu_cores)
{
gpu_cost[1][g] = 0.6;
gpu_cost[2][g] = 0.4;
gpu_cost[3][g] = 0.25;
} // end of for var g scope,

for (var p in cpu_cores)
{
  cpu_cost[1][p] = 0.1;
  cpu_cost[2][p] = 0.08;
  cpu_cost[3][p] = 0.05;
}

  
  for (var m_p in cpu_mem_pages)
{
  cpu_mem_cost[1][m_p] = 0.05;
  cpu_mem_cost[2][m_p] = 0.04;
  cpu_mem_cost[3][m_p] = 0.02;
}

for (var m_g in gpu_mem_pages)
{
  gpu_mem_cost[1][m_g] = 0.15;
  gpu_mem_cost[2][m_g] = 0.12;
  gpu_mem_cost[3][m_g] = 0.06;
}

//////////////////////////

} // End of execute block

///////////////////////////////////////

/// Start to implement the optimization model //////

// Basic Decision variables // 
// Processors 
dvar boolean x1[service_categories][cpu_cores]; // x_{ip}
dvar boolean x2[service_categories][gpu_cores]; // x_{ig}
// memory
dvar boolean x3[service_categories][cpu_mem_pages]; // x_{imp}
dvar boolean x4[service_categories][gpu_mem_pages]; // x_{img}


// Auxilary decision variables 

// Binary variables to gurantee the selection of a cpu or gpu architecture family
dvar boolean y1[service_categories]; 


//////////////////////////////////////////////////////////////////

// Binary decision variables to linearize the multiplication of two decision variables
dvar boolean x1y1[service_categories][cpu_cores];  // x1y1{ip} ... to linearize the multiplication of x_{ip}*y1
dvar boolean x2y1[service_categories][gpu_cores];  // x2y1{ig} ... to linearize the multiplication of x_{ig}*(1-y1) --> x2y1[i][g]
dvar boolean x3y1[service_categories][cpu_mem_pages];  // x3y1{imp} ... to linearize the multiplication of x_{imp}*y1 --> x3y1[i][m_p]
dvar boolean x4y1[service_categories][gpu_mem_pages];  // x4y1{ig} ... to linearize the multiplication of x_{img}*(1-y1) --> x4y1[i][m_g]
/////////////////////////////////

/// Model

minimize ( 
           (sum(i in service_categories) sum(p in cpu_cores)  x1[i][p] * cpu_cost[i][p])  
         + (sum(i in service_categories) sum(g in gpu_cores)  x2[i][g] * gpu_cost[i][g]) 
         + (sum(i in service_categories) sum(m_p in cpu_mem_pages)  x3[i][m_p] * cpu_mem_cost[i][m_p])
         + (sum(i in service_categories) sum(m_g in gpu_mem_pages)  x4[i][m_g] * gpu_mem_cost[i][m_g])
         );

subject to { /// opening constraints
//// starting of processing constraints ////// 
// If f(x) > 0 Then g(x) >= 0
forall (i in service_categories)
	cons027:
	// -g(x) <= M1y1
	(sum(p in cpu_cores) x1[i][p]) <= M1[i]*y1[i];
forall (i in service_categories)
	cons028:
	 // f(x) <= M1(1-y1) 
	 u[i]*ops[i] - delta[i]*v_p[i]*sum(p in cpu_cores) x1y1[i][p] <= M1[i]*(1-y1[i]);
forall (i in service_categories)
	cons0027:
	// -g(x) <= M1y1
	(sum(g in gpu_cores) x2[i][g]) <= M1[i]*(1-y1[i]);
forall (i in service_categories)
	cons0028:	 	
	u[i]*ops[i] - delta[i]*v_g[i]*sum(g in gpu_cores) x2y1[i][g] <= M1[i]*y1[i]; 
///////////////////////////// End of processing constraints ///////////////////////////////////

////////////////////////////////////////////
//// Memory constraints 
////////////////////////////////////////////
forall (i in service_categories)
   consMemCPU:
  	h[i]*n[i]*z - a_p * (sum(m_p in cpu_mem_pages) x3y1[i][m_p]) <= h[i]*n[i]*z*(1-y1[i]);
  	
forall (i in service_categories)
   consMemGPU:
  	h[i]*n[i]*z - a_g * (sum(m_g in gpu_mem_pages) x4y1[i][m_g]) <= h[i]*n[i]*z*y1[i];  		

/////////////////////////////////////////////////////////////////////   
///// End of memory constraints ///////////////////////////
/////////////////////////////////////////////////////////////////////////

////////////////////////////////  Linearization of products ///////////////////////////////// 
/// x1y1
forall (i in service_categories)
	   forall (p in cpu_cores)
	     cons008:
			x1y1[i][p] <= y1[i];
	
	forall (i in service_categories)
	   forall (p in cpu_cores)
	     cons009:
			x1y1[i][p] <= x1[i][p];
	
	forall (i in service_categories)
	   forall (p in cpu_cores)
	     cons010:
			x1y1[i][p] >= x1[i][p] + y1[i] - 1;
			
	forall (i in service_categories)
	   forall (p in cpu_cores)
	     cons011:
			x1y1[i][p] >= 0;

////////////////////////////////////
/// x2y1
forall (i in service_categories)
	   forall (g in gpu_cores)
	     cons8:
			x2y1[i][g] <= 1-y1[i];
	
	forall (i in service_categories)
	   forall (g in gpu_cores)
	     cons9:
			x2y1[i][g] <= x2[i][g];
	
	forall (i in service_categories)
	   forall (g in gpu_cores)
	     cons0010:
			x2y1[i][g] >= x2[i][g] - y1[i];
			
	forall (i in service_categories)
	   forall (g in gpu_cores)
	     cons0011:
			x2y1[i][g] >= 0;
////////////////////////////////////////////////
/// x3y1
forall (i in service_categories)
	   forall (m_p in cpu_mem_pages)
	     consx3008:
			x3y1[i][m_p] <= y1[i];
	
	forall (i in service_categories)
	   forall (m_p in cpu_mem_pages)
	     consx3009:
			x3y1[i][m_p] <= x3[i][m_p];
	
	forall (i in service_categories)
	   forall (m_p in cpu_mem_pages)
	     consx3010:
			x3y1[i][m_p] >= x3[i][m_p] + y1[i] - 1;
			
	forall (i in service_categories)
	   forall (m_p in cpu_mem_pages)
	     consx3011:
			x3y1[i][m_p] >= 0;
////////////////////////////////////
/// x4y1
forall (i in service_categories)
	   forall (m_g in gpu_mem_pages)
	     consx008:
			x4y1[i][m_g] <= 1-y1[i];
	
	forall (i in service_categories)
	   forall (m_g in gpu_mem_pages)
	     consx4009:
			x4y1[i][m_g] <= x4[i][m_g];
	
	forall (i in service_categories)
	   forall (m_g in gpu_mem_pages)
	     cons00010:
			x4y1[i][m_g] >= x4[i][m_g] - y1[i];
			
	forall (i in service_categories)
	   forall (m_g in gpu_mem_pages)
	     cons00011:
			x4y1[i][m_g] >= 0;


///////////////////////////////// End of prodcuts' linearization /////// 
/// isolation constraints ///
forall (p in cpu_cores)
		cons058:
		  (sum(i in service_categories) x1[i][p]) <= 1;  
	forall (g in gpu_cores)
		cons059:
		  (sum(i in service_categories) x2[i][g]) <= 1;
		  
	forall (m_p in cpu_mem_pages) 
		cons060: 
		  (sum(i in service_categories) x3[i][m_p]) <= 1;
		  
	forall (m_g in gpu_mem_pages) 
		cons061:
		  (sum(i in service_categories) x4[i][m_g]) <= 1;
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//// Resource constraints 
	cons062:
	(sum(i in service_categories) sum(p in cpu_cores)  x1[i][p]) <= P;
	cons063:
	(sum(i in service_categories) sum(g in gpu_cores)  x2[i][g]) <= G;
	cons064: /// Relaxed by solver due to M value
	(sum(i in service_categories) sum(m_p in cpu_mem_pages)  x3[i][m_p]) <= M_p;
	cons065: /// Relaxed by solver due to M value
	(sum(i in service_categories) sum(m_g in gpu_mem_pages)  x4[i][m_g]) <= M_g;
} /// closing the constraints 


execute
{  
for (var i in service_categories)
{
  cpu_limit[i] = 0;
  memc_limit[i] = 0;
  gpu_limit[i] = 0;
  memg_limit[i] = 0;
  ///////////////////////////////
  for (var p in cpu_cores)
   {
   	cpu_limit[i] += x1[i][p];  
   }
  for (var m_p in cpu_mem_pages) 
   {
    memc_limit[i] += x3[i][m_p];
   }
  ////////////////////////////////
  
  ///////////////////////////////
  for (var g in gpu_cores)
   {
   	gpu_limit[i] += x2[i][g];  
   }
  for (var m_g in gpu_mem_pages) 
   {
    memg_limit[i] += x4[i][m_g];
   } 
  ///////////////////////////////
 }
}
