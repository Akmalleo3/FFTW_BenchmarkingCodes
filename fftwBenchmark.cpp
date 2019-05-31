#include <utils.h>
#include <fftw3.h>
#include <math.h>
#include <omp.h>
#include <iostream>


/*A Preliminary FFTW benchmark driver that compares
  performance running n fftw executes with either a
  single plan or n plans.  

*/
int main(int argc, char ** argv){

  //Baked in parameters:
  //DOUBLE PRECISION
  //SINGLE DIMENSION - to do

  if(argc < 5){

    std::cout << "Usage: \n ./fftwBenchMark [n_iterations] [FFTW_PLAN_TYPE] [DIMENSIONS] [SIZE (per dimension)]  [n_threads for fftw_execute (default omp_max_threads())] \n for plan types {0:FFTW_MEASURE, 1:FFTW_ESTIMATE, 2: FFTW_PATIENT, 3: FFTW_EXHAUSTIVE}";
    std::cout << std::endl ;
    exit(1);

  }

  unsigned int fftwPlanType;

  int n_iters = std::stoi(argv[1]); 

  std::string planType;
  std::string measure("FFTW_MEASURE");
  std::string estimate("FFTW_ESTIMATE");
  std::string patient("FFTW_PATIENT");
  std::string exhaustive("FFTW_EXHAUSTIVE");

  
  switch(std::stoi(argv[2])){
  case 0:
    fftwPlanType = FFTW_MEASURE;
    planType = measure;
    break;
  case 1:
    fftwPlanType = FFTW_ESTIMATE;
    planType = estimate;
    break;
  case 2:
    fftwPlanType = FFTW_PATIENT;
    planType = patient;
    break;
  case 3:
    fftwPlanType = FFTW_EXHAUSTIVE;
    planType = exhaustive;
    break;
  default: //incorrect usage
    std::cout << "Invalid FFTW Plan type parameter parameter, please try again" << std::endl;
    return 1; 
  }

  int n_dims(std::stoi(argv[3]));
  
  size_t size_per_dim(std::stoi(argv[4]));
  
  int threadsArg {0};
  if(argc == 6){
    threadsArg = std::stoi(argv[5]);
  }

  
  //use number of threads specified, otherwise maximum available
  int max_threads =  omp_get_max_threads();
  int n_threads = threadsArg == 0 ? max_threads : threadsArg;
  if (n_threads > max_threads)
    n_threads = max_threads; //ensure user did not specify too many threads
    
  
  //unless user specified single threaded
  if(n_threads != 1){ 
  //prepare for multithreaded fftw
  int err = fftw_init_threads();
    if(!err){
      std::cout << "fftw_init_threads failed" << std::endl;
      return 1;
    }
    
    fftw_plan_with_nthreads(n_threads);
    std::cout << "calls to fftw_execute will utilize " << n_threads << " threads" << std::endl;
  }
  
  FFTW_FORGET_WISDOM();

  CNTime timer; 

  fftw_complex *in, *out;
  fftw_plan p;

  size_t contiguousArraySize;
  int n[3]{0,0,0}; 
  switch(n_dims){
  case 1:
    contiguousArraySize = size_per_dim;
    n[0] = size_per_dim;
    break;
  case 2:
    contiguousArraySize = size_per_dim*size_per_dim;
    n[0] = size_per_dim;
    n[1] = size_per_dim;
    break;
  case 3:
    contiguousArraySize = size_per_dim*size_per_dim*size_per_dim;
    n[0] = size_per_dim;
    n[1] = size_per_dim;
    n[2] = size_per_dim;
    break;
  default:
    std::cout << "Invalid dimension parameter, please try again" << std::endl;
    return 1; 
  }
  
  in = fftw_alloc_complex(contiguousArraySize); 
  out = fftw_alloc_complex(contiguousArraySize);

  std::cout << "*---------------------------------------------------------------------*" << std::endl;
  std::cout << "Executing " << n_iters << " iterations with a shared plan of type " << planType << std::endl; 


  timer.start();
  p = fftw_plan_dft(n_dims, n, in, out, FFTW_FORWARD, fftwPlanType);
  double planConstruction1 = timer.elapsedsec();
    
  std::cout << "Time for plan construction: " << planConstruction1 << std::endl;

  unsigned int dummy = 1;
  //initialize data after initializing the plan
  for(int i = 0; i < contiguousArraySize; i++){
    in[i][0] = M_PI*randm11r(&dummy);
    in[i][1] = M_PI*randm11r(&dummy);
  }
  

  double executeTime;
  double totalExecTime1{0};
  for(int i = 0; i < n_iters; i++){
    timer.restart();
    fftw_execute(p);
    executeTime = timer.elapsedsec();
    totalExecTime1 += executeTime;

  }
  std::cout << "Total Execution Time: " << totalExecTime1 << std::endl; 
  std::cout << "Final Total for this trial: " << planConstruction1 + totalExecTime1 << std::endl;
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out); // no error checking yet



  //Remove the advantage of collected wisdom
  FFTW_FORGET_WISDOM();

  std::cout << "*---------------------------------------------------------------------*" << std::endl;
  std::cout << "Executing " << n_iters << " iterations, initializing a new " << planType << " plan each time " << std::endl;


  if(threadsArg != 1)
    fftw_plan_with_nthreads(n_threads);

  fftw_complex *in2, *out2;
  fftw_plan p2;
  in2 = fftw_alloc_complex(contiguousArraySize); 
  out2 = fftw_alloc_complex(contiguousArraySize);

  double planConstruction2;
  double totalPlanConstruction2{0};
  double totalExecTime2{0};

  for(int i = 0; i < n_iters; i++){
    timer.restart();
    p2 = fftw_plan_dft(n_dims, n, in2, out2, FFTW_FORWARD, fftwPlanType);
    planConstruction2 = timer.elapsedsec();
    if(i == 0){
      std::cout << "Time in first call to construct plan: " << planConstruction2 << std::endl;
     }
    totalPlanConstruction2 += planConstruction2; 

    //initialize data for this round
    for(int k = 0; k < contiguousArraySize; k++){
      in2[k][0] = M_PI*randm11r(&dummy);
      in2[k][1] = M_PI*randm11r(&dummy);
    }

    timer.restart();
    fftw_execute(p2);
    executeTime = timer.elapsedsec();
    totalExecTime2 += executeTime; 
    
  }
  
  std::cout << "Total Time for Construction of all plans: " << totalPlanConstruction2 << std::endl;
  std::cout << "Total Time for Execution " << totalExecTime2 << std::endl;
  std::cout << "Final Total for this Trial: " << totalPlanConstruction2 + totalExecTime2 << std::endl;
  
  fftw_destroy_plan(p2);
  fftw_free(in2); fftw_free(out2);
  fftw_cleanup_threads();
  FFTW_FORGET_WISDOM();

}

