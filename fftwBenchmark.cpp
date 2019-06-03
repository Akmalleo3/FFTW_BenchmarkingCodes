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

  //allocate result array and set up dimension array
  size_t arraySize;
  int dims[3]{0,0,0}; 
  int idist; // parameter for plan_many
  int odist; //parameter for plan_many

  switch(n_dims){
  case 1:
    arraySize = size_per_dim;
    dims[0] = size_per_dim;
    idist = dims[0];
    odist = dims[0];
    break;
  case 2:
    arraySize = size_per_dim*size_per_dim;
    dims[0] = size_per_dim;
    dims[1] = size_per_dim;
    idist = dims[0]*dims[1]; 
    odist = dims[0]*dims[1]; 
    break;
  case 3:
    arraySize = size_per_dim*size_per_dim*size_per_dim;
    dims[0] = size_per_dim;
    dims[1] = size_per_dim;
    dims[2] = size_per_dim;
    idist = dims[0]*dims[1]*dims[2];
    odist = dims[0]*dims[1]*dims[2];
    break;
  default:
    std::cout << "Invalid dimension parameter, please try again" << std::endl;
    return 1; 
  }
  
  in = fftw_alloc_complex(arraySize); 
  out = fftw_alloc_complex(arraySize);

  std::cout << "*---------------------------------------------------------------------*" << std::endl;
  std::cout << "Executing " << n_iters << " iterations with a shared plan of type " << planType << std::endl; 


  timer.start();
  p = fftw_plan_dft(n_dims, dims, in, out, FFTW_FORWARD, fftwPlanType);
  double planConstruction1 = timer.elapsedsec();
    
  std::cout << "Time for plan construction: " << planConstruction1 << std::endl;

  unsigned int dummy = 1;
  //initialize data after initializing the plan
  for(int i = 0; i < arraySize; i++){
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
  double total1 = planConstruction1 + totalExecTime1;
  std::cout << "Total for this expt 1: " << total1 << std::endl;
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
  in2 = fftw_alloc_complex(arraySize); 
  out2 = fftw_alloc_complex(arraySize);

  double planConstruction2;
  double planRepeated2=0.0, regendata2 = 0.0;
  double totalExecTime2{0};

  for(int i = 0; i < n_iters; i++){
    timer.restart();
    p2 = fftw_plan_dft(n_dims, dims, in2, out2, FFTW_FORWARD, fftwPlanType);
    if(i == 0) {
      planConstruction2 = timer.elapsedsec();
    } else
      planRepeated2 += timer.elapsedsec(); 

    //initialize data for this round
    timer.restart();
    for(int k = 0; k < arraySize; k++){
      in2[k][0] = M_PI*randm11r(&dummy);
      in2[k][1] = M_PI*randm11r(&dummy);
    }
    regendata2 += timer.elapsedsec(); 

    timer.restart();
    fftw_execute(p2);
    executeTime = timer.elapsedsec();
    totalExecTime2 += executeTime; 
    
  }
  
  std::cout << "First plan Construction: " << planConstruction2 << std::endl;
  std::cout << "Total for all later plans: " << planRepeated2 << std::endl;
  std::cout << "Mean time for later plan (us): " << 1e6*planRepeated2/(n_iters-1) << std::endl;
  std::cout << "Total Time for Execution " << totalExecTime2 << std::endl;
  double total2 = planRepeated2+totalExecTime2+planConstruction2+regendata2;
  std::cout << "Total time for expt 2: " << total2 << std::endl;

  
  fftw_destroy_plan(p2);
  fftw_free(in2); fftw_free(out2);
  
  FFTW_FORGET_WISDOM();


  int totalSize = arraySize*n_iters;
  
  std::cout << "*-------------------------------------------------------------------------*" << std::endl;
  std::cout << "Executing " << n_iters << " iterations, sharing a single fftw_plan_many_dft of type" << planType  << std::endl;
  
  int istride = 1;
  int ostride = 1; //array is contiguous in memory


  std::cout << "Attempting to allocate an array of " << totalSize << "  elements" << std::endl; 

  //allocate space for in and out arrays for all n_iter executions
  //contiguous in memory
  fftw_complex *in3, *out3;
  fftw_plan p3;
  in3 = fftw_alloc_complex(arraySize*n_iters); 
  out3 = fftw_alloc_complex(arraySize*n_iters);

  if(in3 == NULL | out3 == NULL){

    std::cout << "program parameters request too large a memory allocation for trial 3.  \n Aim for fewer iterations and smaller points in each dimension" << std::endl; 
  return 1;
  }

  std::cout << "SUCCESS" << std::endl;
  
  timer.start();
  p3 = fftw_plan_many_dft(n_dims, dims, n_iters, in3, dims, istride, idist, out3, dims, ostride, odist, FFTW_FORWARD, fftwPlanType);
  double planConstruction3 = timer.elapsedsec();
    
  std::cout << "Time for many_plan construction: " << planConstruction3 << std::endl;

    //initialize data after initializing the plan
  for(int i = 0; i < totalSize; i++){
    in3[i][0] = M_PI*randm11r(&dummy);
    in3[i][1] = M_PI*randm11r(&dummy);
  }

  timer.restart();
  fftw_execute(p3);
  double totalExecTime3 = timer.elapsedsec();

  std::cout << "Total Execution Time: " << totalExecTime3 << std::endl;
  double total3 = planConstruction3 + totalExecTime3;
  std::cout << "Total time for expt 3: " << total3 << std::endl;

  std::cout << "*----------------------------------------------------------*" << std::endl;
  std::cout << "Total time : " << total1+total2+total3 << std::endl;

  fftw_destroy_plan(p3);
  fftw_free(in3); fftw_free(out3); // no error checking yet


  FFTW_FORGET_WISDOM();
  fftw_cleanup_threads();
}

