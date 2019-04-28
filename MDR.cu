#include <iostream>
#include <math.h>
#include <time.h>
//#include "MDR.h"
//#include "MDR_kernel.cu"
//#include "MDR.cu"

#if _WIN32
    //Windows threads.
    #include <windows.h>

    typedef HANDLE CUTThread;
    typedef unsigned (WINAPI *CUT_THREADROUTINE)(void *);

    #define CUT_THREADPROC unsigned WINAPI
    #define  CUT_THREADEND return 0

#else
    //POSIX threads.
    #include <pthread.h>

    typedef pthread_t CUTThread;
    typedef void *(*CUT_THREADROUTINE)(void *);

    #define CUT_THREADPROC void
    #define  CUT_THREADEND
#endif

//Create thread.
CUTThread start_thread( CUT_THREADROUTINE, void *data );

//Wait for thread to finish.
void end_thread( CUTThread thread );

//Destroy thread.
void destroy_thread( CUTThread thread );

//Wait for multiple threads.
void wait_for_threads( const CUTThread *threads, int num );

#if _WIN32
    //Create thread
    CUTThread start_thread(CUT_THREADROUTINE func, void *data){
        return CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, data, 0, NULL);
    }

    //Wait for thread to finish
    void end_thread(CUTThread thread){
        WaitForSingleObject(thread, INFINITE);
        CloseHandle(thread);
    }

    //Destroy thread
    void destroy_thread( CUTThread thread ){
        TerminateThread(thread, 0);
        CloseHandle(thread);
    }

    //Wait for multiple threads
    void wait_for_threads(const CUTThread * threads, int num){
        WaitForMultipleObjects(num, threads, true, INFINITE);

        for(int i = 0; i < num; i++)
            CloseHandle(threads[i]);
    }

#else
    //Create thread
    CUTThread start_thread(CUT_THREADROUTINE func, void * data){
        pthread_t thread;
        pthread_create(&thread, NULL, func, data);
        return thread;
    }

    //Wait for thread to finish
    void end_thread(CUTThread thread){
        pthread_join(thread, NULL);
    }

    //Destroy thread
    void destroy_thread( CUTThread thread ){
        pthread_cancel(thread);
    }

    //Wait for multiple threads
    void wait_for_threads(const CUTThread * threads, int num){
        for(int i = 0; i < num; i++)
            end_thread( threads[i] );
    }

#endif




#define imin(a,b) (a<b?a:b)
#define imax(a,b) (a>b?a:b)
#define func(a)

char* phenoFile;
char* genoFile;
char* outputFile;
char* combFile;
int basic_model;

float THR = THR ; //leave this space
#define NSNPS NSNPS
#define NUMCOMBS NUMCOMBS
#define CUT CUT
#define NIND NIND
#define BSx BSx
#define ORDER ORDER
#define CV CV
#define TABLE_SIZE TABLE_SIZE
#define NUMDEVICES NUMDEVICES
#define ONEORTWO ONEORTWO //two if CV > 1
#define GSx ((NUMCOMBS/NUMDEVICES+BSx-1) / BSx )
#define MEASURE 'MEASURE'

#define mat_SNP_size NIND * NSNPS * sizeof(int)
#define v_pheno_size NIND * sizeof(int)
#define output_size NUMCOMBS *  ONEORTWO * CV * sizeof(float) //oneotwo: 2 is one for train and one for test
#define fp_size NUMCOMBS * sizeof(int)
#define tp_size NUMCOMBS * sizeof(int)
#define combinations_size NUMCOMBS * ORDER * sizeof(int)
#define indices_size NIND * sizeof(int)


#define TESTCOMB -112
#define TESTSNP0 -87509
#define TESTSNP1 -370675

struct controlscases {
int controls;
int cases;
};


 struct str
{
   float value;
   int index;
 };

struct DataStruct {
    int		deviceID;
    int		deviceCount;
    int*	mat_SNP;
    int*	combinations;
    int* 	v_pheno;
    int* 	cv_indices;
    float*	output;
    int*	tp;
    int*	fp;
    float	start_clock;
};


__device__ float compute_measure(int cases_high, int controls_high, int controls_low, int cases_low, char m){
	float train_measure;
	//int positives = cases_high + cases_low
	//int negatives = controls_low + controls_high
	if (m - '0' == 50){ //BALANCED ACCURACY
		train_measure = (cases_high/float(cases_high + cases_low) + controls_low/float(controls_low + controls_high))/2;
	}
	else if (m - '0' == 49){ //ACCURACY
		train_measure = float(cases_high + controls_low)/float(cases_high + controls_low + controls_high + cases_low);
		
	}
	else if (m - '0' == 55){ //TODO
		train_measure = float(cases_high + controls_low)/float(cases_high + controls_low + controls_high + cases_low);
		
	}
	return train_measure;
}

__device__ int dev_pow(int b,int e) {
	int o = 1;
	if (e == 0)
		return 1;
	for (int i=0; i<e; i++){
		o *= b;
	}
	return o;
}

__device__ int count_digits(int i) {
	if (i < 10)
		return 1;
	if (i < 100)
		return 2;
	if (i < 1000)
		return 3;
	if (i < 10000)
		return 4;
	if (i < 100000)
		return 5;
	return 6;
}

//base decimal to base 3
__device__ void int_to_index(int n, int order, int* v){
	int r;
	for(int i=0; i<order; i++){
		float num = (float)n;
		int den = dev_pow(3,order-i-1);
		r = (int)(num/den);
		n -= r * dev_pow(3,order-i-1);
		v[i] = r ;
	}
	return;
}

//base 3 to base decimal
__device__ int index_to_int(int* v_indices, int order){
	int o = 0;
	for(int i=0; i<order; i++){
		o += v_indices[i] * dev_pow(3, order-i-1);
	}
	return o;
}


//#include "MDR.h"
//#include "MDR.cu"

__constant__ int dev_v_pheno[NIND];
__constant__ int dev_cv_indices[NIND];
	

__global__ void MDR( int* dev_SNP_values, float* dev_output, int* dev_tp, int* dev_fp, int* dev_combinations, float THR, int deviceID, int deviceCount) {
    
    	
	//printf(" %d + %d * %d :", threadIdx.x, blockIdx.x, blockDim.x);
	//__shared__ float cache[BS][threadsPerBlock];
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	
	
	
	//printf(" %d ", tid);
	int d;
	if (deviceID + 1 < deviceCount)
		d = (NUMCOMBS/deviceCount);
	else
		d = NUMCOMBS - 1 - (((0 + ((NUMCOMBS *  ONEORTWO * CV)/(deviceCount) * deviceID ))) - 1); //how many total - how many done
	
	if (tid >= d)
		return;

	//printf(" %d ", tid);
	//int* thread_combination = (int*)malloc(ORDER * sizeof(int));
	int thread_combination[ORDER]; //a combination (thread level)
	//retrieve the combination indices
	for (int i=0; i< ORDER; i++) {
		*(&thread_combination[0] + i) = *(dev_combinations + tid * ORDER + i);
	}
	
	
	//printf("thread with tid %d is assigned combination: <%d, %d>\n", tid, thread_combination[0], thread_combination[1]); 
	
	//retrieve the genotype of each snp in the combination, from SNPvalues, for ALL individuals
	//int* thread_geno = (int*)malloc(NIND * ORDER * sizeof(int));
	/*
	int thread_geno[ORDER * NIND];
	for (int i=0; i< NIND; i++) {
		for (int j=0; j< ORDER; j++) {
			*(&thread_geno[0] + j * NIND + i  ) = *(dev_SNP_values + NIND * *(&thread_combination[0] + j) + i);
		}
	}
	*/
	
	/*
	if (*(&thread_combination[0] + 0) == TESTSNP0 && *(&thread_combination[0] + 1) == TESTSNP1){
		printf("START of SNP_values\n");
		for (int j=0; j < 4050; j++) {
			printf("%d ", *(dev_SNP_values + NIND * 0 + j));
			if (j == 1999 || j == 4000)
				printf("\n");
			}
		printf("\nend START\n\n");
		}
		
	*/
	
	/*
	if (*(&thread_combination[0] + 0) == TESTSNP0 && *(&thread_combination[0] + 1) == TESTSNP1){
		printf("GENOs for %d and %d\n", TESTSNP0, TESTSNP1);
		for (int j=0; j < ORDER; j++) {
			printf("\n");
			for (int i=0; i < NIND; i++) {
				if (i< 10 || i > NIND-10)
					printf("%d ", *(dev_SNP_values + NIND * *(&thread_combination[0] + j) + i));
			}
		}
		printf("\nend GENOs");
	}
	*/
	
	
	struct controlscases thread_table[TABLE_SIZE];

	//replace this initialization?
	for (int i=0; i< TABLE_SIZE; i++) {
		(*(&thread_table[0] + i )).controls = 0;
		(*(&thread_table[0] + i )).cases = 0;
		}
		
	int geno[ORDER]; //support variable for counting, stores a geno combination
	int index_in_table = 0; //from a geno, to the index in the table
	int indiv;
	int cases_high;
	int controls_high;
	int cases_low;
	int controls_low;
	// ba = ((tp/p) + (tn/n))/2 = (sensitivity + specificity) /2
	float train_measure;
	float test_measure;
	int v[ORDER];
	int ph;
	//CV loop
	for (int cv=0; cv<CV; cv++){

		
		
		//*****************
		//TRAINING
		//*****************
		//populate the 3^ORDER-tot-entries table
		for (int n=0; n< NIND; n++) {
			if (CV > 1){
				if ((n >= int((cv/float(CV))*NIND)) && (n <= int(((cv+1)/float(CV))*NIND)) )//reserved for testing
				 		continue;
				 }
			 indiv = *(&dev_cv_indices[0] + n);
			 for (int i=0; i< ORDER; i++) 
			 	geno[i] = *(dev_SNP_values + NIND * *(&thread_combination[0] + i) + indiv); //i-th snp geno
			 index_in_table = index_to_int(geno, ORDER);
			 
			 if (int(*(dev_v_pheno + indiv))) { //get the pheno
			 	
			 	
			 	(*(&thread_table[0] + index_in_table )).cases += 1;
			 }
			 else{
			 	//if (tid == TESTCOMB)
			 	//	printf("?????? geno %d and healthy ph: %d\n",index_in_table, indiv);
			 	(*(&thread_table[0] + index_in_table )).controls += 1;
			 }
		}

		//only a print
		if (tid == TESTCOMB || (*(&thread_combination[0] + 0) == TESTSNP0 && *(&thread_combination[0] + 1) == TESTSNP1)){
			printf("***************\ngpu%d-tid%d\ncomb. ", deviceID, tid);
			for (int q=0; q< ORDER; q++)
				printf("%d ", *(&thread_combination[0] + q));
			printf("\n\n");
			for (int i=0; i< TABLE_SIZE; i++) {
				printf("thread_table[%d].controls, cases: %d %d ",i,(*(&thread_table[0] + i )).controls, (*(&thread_table[0] + i )).cases);
				if ( (((*(&thread_table[0] + i )).cases) / float((*(&thread_table[0] + i )).controls + 0.01) >= THR ) ){
					int_to_index(i, ORDER, v);
					printf(" geno ");
					for (int l=0; l< ORDER; l++)
						printf("%d ", v[l] );
					printf("is HIGH\n");
				}
				else
					printf("\n");
			
			}
			printf("\n");
		}
	
		//moving two a two-dim variable
		cases_high = 0;
		controls_high = 0;
		cases_low = 0;
		controls_low = 0;
		//int c = 0;
		for (int i=0; i< TABLE_SIZE; i++) {
			int_to_index(i, ORDER, v);
				
			if ( (((*(&thread_table[0] + i )).cases) / float((*(&thread_table[0] + i )).controls + 0.01) >= THR )) {
				cases_high += (*(&thread_table[0] + i )).cases;
				controls_high += (*(&thread_table[0] + i )).controls;

			}
			else{
				//here in LOW also the case 0 controls 0 cases
				cases_low += (*(&thread_table[0] + i )).cases;
				controls_low += (*(&thread_table[0] + i )).controls;
				if (tid == TESTCOMB){
					printf("tid %d (comb. ", tid) ;
					for (int q=0; q< ORDER; q++)
						printf("%d ", *(&thread_combination[0] + q));
					printf("), geno ");
					for (int l=0; l< ORDER; l++)
						printf("%d ", v[l] );
					printf("is LOW\n");
				}

			}
		}
		/*
		if (CV > 1)
			*(&high_genos[0] + c*ORDER + 0) = 9; //end sequence, since high_genos only reports the high ones
		*/
		//printf("******************\n");
		
		train_measure = compute_measure(cases_high, controls_high, controls_low, cases_low, MEASURE);
		
		
		
		//only a print
		if (tid == TESTCOMB || (*(&thread_combination[0] + 0) == TESTSNP0 && *(&thread_combination[0] + 1) == TESTSNP1)){
		if (MEASURE - '0' == 50){
		printf("(tid %d) TRAIN BA %1.5f = (%d/(%d+%d) + %d/(%d+%d))/2\n***************\n", 
				 tid, train_measure, cases_high, cases_high, controls_high, controls_low, controls_low, cases_low);
		}
		else if (MEASURE - '0' == 49){
		printf("(tid %d) TRAIN AC	 %1.5f = (%d+%d)/(%d+%d+%d+%d)\n***************\n", 
			 tid, train_measure, cases_high, controls_low, cases_high, controls_low, controls_high, cases_low);
		}
		else if (MEASURE - '0' == 55){
		printf("(tid %d) TRAIN AC	 %1.5f = (%d+%d)/(%d+%d+%d+%d)\n***************\n", 
			 tid, train_measure, cases_high, controls_low, cases_high, controls_low, controls_high, cases_low);
		
		}
		}
		
		
		
		
		

		//write result to global memory
		if (CV == 1){
			*(dev_output + NUMCOMBS * 0 + 1 * tid + 0) = train_measure;
			*(dev_tp + NUMCOMBS * 0 + 1 * tid + 0) = cases_high;
			*(dev_fp + NUMCOMBS * 0 + 1 * tid + 0) = controls_high;
			}
		else
			*(dev_output + NUMCOMBS * cv + 2 * tid + 0) = train_measure;
			//TODO
		
		
		if (CV > 1)
		{ //TODO
			//*****************
			//TESTING
			//*****************
		
			if (tid == TESTCOMB) printf("CV-TEST %d/%d:\n", cv+1, CV);
		
			cases_high = 0;
			controls_high = 0;
			cases_low = 0;
			controls_low = 0;
			for (int n=0; n< NIND; n++) {
				 if ((n < int((cv/float(CV))*NIND)) || (n > int(((cv+1)/float(CV))*NIND)) )//reserved for training
				 	continue;
				 
				 indiv = *(&dev_cv_indices[0] + n);
				 for (int i=0; i< ORDER; i++) 
				 	geno[i] = *(dev_SNP_values + NIND * *(&thread_combination[0] + i) + indiv); //i-th snp geno
				 ph = int(*(dev_v_pheno + indiv));
				 
				 //check if retrieved geno is in high or low
				 index_in_table = index_to_int(geno, ORDER);
				 if (((*(&thread_table[0] + index_in_table )).cases)/((*(&thread_table[0] + index_in_table )).controls + 0.01) >= THR){
				 	if (ph)
						cases_high += 1;
					else
						controls_high += 1;
				 }
				 else{
				 	if (ph)
						cases_low += 1;
					else
						controls_low += 1;
				 }
				 	
				 
				 
				 /*
				 for (int i=0; i< TABLE_SIZE * ORDER; i++) {
				 	 if (high_genos[i] == 9){ //reached the end
				 	 	if (ph)
						 	cases_low += 1;
						 else
						 	controls_low += 1;
						 break;
						 }
				 	 	
				 	 int isequal = 1;
				 	 for (int j=0; j< ORDER; j++){
					 	if  (high_genos[i + j] != geno[i])
					 		isequal = 0;
					 		break;
					 	}
					 if (isequal){
					 	if (ph)
					 		cases_high += 1;
					 	else
					 		controls_high += 1;
					 	break; //found, exit loop;
					 	
					 }
			
					 	 
				 }
				 */
			}
	
		
		

			//test_measure = float(controls_high + cases_low)/float(cases_high + controls_high + cases_low + controls_low);
			if (MEASURE - '0' == 50){
			if (controls_low + cases_low == 0)
				//train_measure = float(controls_high + cases_low)/float(cases_high + controls_high + cases_low + controls_low);
				test_measure = (cases_high/float(cases_high + controls_high) + 0)/1;
			
			else if (cases_high + controls_high == 0)
				test_measure = (0 + controls_low/float(controls_low + cases_low))/1;
			else
				test_measure = (cases_high/float(cases_high + controls_high) + controls_low/float(controls_low + cases_low))/2;
		
			if (tid == TESTCOMB || (*(&thread_combination[0] + 0) == TESTSNP0 && *(&thread_combination[0] + 1) == TESTSNP1))
			printf("(tid %d) TRAIN BA %1.5f = (%d/(%d+%d) + %d/(%d+%d))/2\n", 
					 tid, test_measure, cases_high, cases_high, controls_high, controls_low, controls_low, cases_low);
			}
			else if (MEASURE - '0' == 49){
				test_measure = float(cases_high + controls_low)/float(cases_high + controls_low + controls_high + cases_low);
				if (tid == TESTCOMB || (*(&thread_combination[0] + 0) == TESTSNP0 && *(&thread_combination[0] + 1) == TESTSNP1))
				printf("(tid %d) TRAIN AC	 %1.5f = (%d+%d)/(%d+%d+%d+%d)\n", 
					 tid, test_measure, cases_high, controls_low, cases_high, controls_low, controls_high, cases_low);
			
			}
			//write result to global memory
			*(dev_output + NUMCOMBS * cv + 2 * tid + 1) = test_measure;
		
			 
			if (tid == TESTCOMB) printf("**********************************\n\n"); 
		
		}
	}
	
}

float* extract_min(struct str** top_cut_list, int cut, float* a){
	float minim = 1;
	float ind = -1;
	for (int i = 0; i < cut; i++){
		//fprintf(stderr,"inside func: value,index: %f : %d\n", (*(top_cut_list+i))-> value, (*(top_cut_list+i))-> index);
		if ((*(top_cut_list+i))-> value < minim){
			minim = (*(top_cut_list+i))-> value;
			ind = i;
		}
	}
	a[0] = minim;
	a[1] = ind;
	//fprintf(stderr,"%f : %f\n", *(a+0), *(a+1));
	return a;
}

void merge(struct str** arr, int l, int m, int r) 
{ 
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 =  r - m; 
  
    /* create temp arrays */
    struct str* L = (struct str*) malloc(n1 * sizeof(struct str));
    struct str* R = (struct str*) malloc(n2 * sizeof(struct str)); 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) 
        *(L + i) = *(*(arr + l + i)); 
    for (j = 0; j < n2; j++) 
        *(R + j) = *(*(arr + m + 1+ j));
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2) 
    { 
        if ((L+i)->value > (R + j)->value) 
        { 
            *(*(arr + k)) = *(L+i); 
            i++; 
        } 
        else
        { 
            *(*(arr + k)) = *(R+j); 
            j++; 
        } 
        k++; 
    } 
  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 
        *(*(arr + k)) = *(L+i);
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 
        *(*(arr + k)) = *(R+j);  
        j++; 
        k++; 
    } 
    free(L);
    free(R);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort(struct str** arr, int l, int r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        int m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort(arr, l, m); 
        mergeSort(arr, m+1, r); 
  
        merge(arr, l, m, r); 
    } 
} 



void parseArgs(int argc, char **argv){
  int i=1;
  if(argc <= 3){
  
    printf("\nmissing input!\n");
	printf("\n\n");
    exit(0);
  }

  while(i<argc){
    if(!strcmp(argv[i], "-gf"))
      genoFile = argv[++i];
    else if(!strcmp(argv[i], "-cf"))
      combFile = argv[++i];
    else if(!strcmp(argv[i], "-pf"))
      phenoFile = argv[++i];
    else if(!strcmp(argv[i], "-out"))
      outputFile = argv[++i];
    else if(!strcmp(argv[i], "-basic_model"))
      basic_model = atoi(argv[++i]);
    else{
      printf("%s : argument not valid! \n",argv[i]);
      exit(1);
    }
    i++;
  }

  if( !genoFile || !phenoFile || !combFile || !outputFile){
    printf("no files specified	.. exiting\n");
    exit(1);
  }
  return;

}


static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "ERROR HANDLED: %s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
  	}
  }
  
void readintData(char *dataFile, unsigned int rows, unsigned int cols, int * data){
  FILE *fp;
  int *dp = data;
  int i;

  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
    //exit(1);
  } 
  
  for (i=0; i<rows*cols; ++i){
	  fscanf(fp, "%d", dp);
	  dp++;
  } 
  fclose(fp);
  return;
}


void readCombinations(char *dataFile, int rows, int cols, int * data){
  FILE *fp;
  int *dp = data;
  int i;

  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
    //exit(1);
  } 
  
  for (i=0; i<rows*cols; ++i){
	  fscanf(fp, "%d", dp);
	  dp++;
  } 
  fclose(fp);
  return;
}
  	
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


#define HANDLE_NULL( a ) {if (a == NULL) { printf( "Host memory failed in %s at line %d\n",  __FILE__, __LINE__ ); ( EXIT_FAILURE );}}

void print_cudaGetDeviceProperties(){
	cudaDeviceProp prop;
	int count;
	HANDLE_ERROR( cudaGetDeviceCount( &count ) );
	
	if (count == 0) {
        	fprintf(stderr,"error in print_cudaGetDeviceProperties: no devices supporting CUDA.\n");
        	return;
    	}
	
	for (int i=0; i< count; i++) {
		HANDLE_ERROR( cudaGetDeviceProperties( &prop, i ) );
		printf( "   --- General Information for device %d ---\n", i );
		printf( "Name:  %s\n", prop.name );
		printf( "Compute capability:  %d.%d\n", prop.major, prop.minor );
		printf( "Clock rate:  %d\n", prop.clockRate );
		printf( "Device copy overlap:  " );
		if (prop.deviceOverlap)
		    printf( "Enabled\n" );
		else
		    printf( "Disabled\n");
		printf( "Kernel execution timeout :  " );
		if (prop.kernelExecTimeoutEnabled)
		    printf( "Enabled\n" );
		else
		    printf( "Disabled\n" );

		printf( "   --- Memory Information for device %d ---\n", i );
		printf( "Total global mem:  %ld\n", prop.totalGlobalMem );
		printf( "Total constant Mem:  %ld\n", prop.totalConstMem );
		printf( "Max mem pitch:  %ld\n", prop.memPitch );
		printf( "Texture Alignment:  %ld\n", prop.textureAlignment );

		printf( "   --- MP Information for device %d ---\n", i );
		printf( "Multiprocessor count:  %d\n",
		            prop.multiProcessorCount );
		printf( "Shared mem per mp:  %ld\n", prop.sharedMemPerBlock );
		printf( "Registers per mp:  %d\n", prop.regsPerBlock );
		printf( "Threads in warp:  %d\n", prop.warpSize );
		printf( "Max threads per block:  %d\n",
		            prop.maxThreadsPerBlock );
		printf( "Max thread dimensions:  (%d, %d, %d)\n",
		            prop.maxThreadsDim[0], prop.maxThreadsDim[1],
		            prop.maxThreadsDim[2] );
		printf( "Max grid dimensions:  (%d, %d, %d)\n",
		            prop.maxGridSize[0], prop.maxGridSize[1],
		            prop.maxGridSize[2] );
		printf( "\n" );
	}
	return;
}



void* routine( void *pvoidData) {
	DataStruct *data = (DataStruct*)pvoidData;
	HANDLE_ERROR( cudaSetDevice( data->deviceID ) );
	
	int deviceID = data->deviceID;
	int deviceCount = data->deviceCount;
  	int* dev_mat_SNP;
  	int* dev_combinations;
	float* dev_output;
	int* dev_fp;
	int* dev_tp;
	
	int* mat_SNP = data->mat_SNP;
	int* combinations = data->combinations;
	int* v_pheno = data->v_pheno;
	int* cv_indices = data->cv_indices;
	float* output = data->output;
	int* tp = data->tp;
	int* fp = data->fp;
	clock_t start_clock = data->start_clock;
 	
  	//Allocate device memory
  	cudaMalloc((void**)&dev_mat_SNP, mat_SNP_size);
  	
	//cudaMalloc((void**)&dev_v_pheno, dev_v_pheno.mem_size); //no need, constant mem!
  	//cudaMalloc((void**)&dev_cv_indices, indices_size); //no need, constant mem!
  	
  	// Copy host memory to device
  	//HANDLE_ERROR( cudaMemcpy(dev_v_pheno, v_pheno, dev_v_pheno.mem_size, cudaMemcpyHostToDevice));
  	cudaMemcpyToSymbol( dev_v_pheno,  v_pheno,  v_pheno_size );
  	cudaMemcpyToSymbol( dev_cv_indices,  cv_indices,  indices_size );
	HANDLE_ERROR( cudaMemcpy(dev_mat_SNP, mat_SNP, mat_SNP_size, cudaMemcpyHostToDevice));
	// (combinations_size / sizeof(int))/deviceCount) * deviceID is the same for all gpus, except last one is more loaded
	
	
	//start index (included), number of elements , end index (included)
	unsigned long s,d, e;
	s = (0 + (NUMCOMBS/deviceCount) * deviceID );
	if (deviceID + 1 < deviceCount)
		d = (NUMCOMBS/deviceCount);
	else
		d = NUMCOMBS - 1 - (s - 1); //how many total - how many done
	e = s + d - 1;
	cudaMalloc((void**)&dev_combinations, (d *  ORDER * sizeof(int)) );
	HANDLE_ERROR( cudaMemcpy(dev_combinations, (combinations + s * ORDER), (d *  ORDER * sizeof(int)), cudaMemcpyHostToDevice));
	
	fprintf(stderr,"\nGPU %d, calling the kernel with this configuration:\n", deviceID);
	fprintf(stderr," comb-start: %lu, #combs: %lu, comb-end: %lu\n order: %d\n NSNPS: %d\n NIND: %d\n # CVs: %d\n THRESHOLD: %f\n BLOCK SIZE: %d\n GRID SIZE: %d\n",s,d,e, ORDER, NSNPS, NIND, CV, THR, BSx, GSx);

  	//HANDLE_ERROR( cudaMemcpy(dev_cv_indices, cv_indices, indices_size, cudaMemcpyHostToDevice));
  	//fprintf(stderr,"matrices copied  to GPU %d\n", deviceID);
  	
  	//cudaHostAlloc((void**)&output,output.mem_size,cudaHostAllocDefault);
  	
  	
  	s = (0 + ((NUMCOMBS *  ONEORTWO * CV)/(deviceCount) * deviceID ));
	if (deviceID + 1 < deviceCount)
		d = (NUMCOMBS *  ONEORTWO * CV)/(deviceCount);
	else
		d = (NUMCOMBS *  ONEORTWO * CV) - 1 - (s - 1); //how many total - how many done
	e = s + d - 1;
	//fprintf(stderr,"GPU %d, out-start: %lu, #out-values: %lu, out-end: %lu\n",deviceID,s,d,e);
	
	cudaMalloc((void**)&dev_output, d * sizeof(float));
  	cudaMalloc((void**)&dev_tp, d * sizeof(int));
  	cudaMalloc((void**)&dev_fp, d * sizeof(int));
  	
	// kernel call
	dim3 dimBlock(BSx);//,BSy,BSz);
	dim3 dimGrid(GSx);//,GSy,GSz);
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float elapsedTime;
	cudaEventRecord(start, 0);
	
	//fprintf(stderr,"GPU %d, calling kernel; elapsed %f seconds\n", deviceID, ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
	
	MDR<<< dimGrid, dimBlock >>>(dev_mat_SNP, dev_output, dev_tp, dev_fp, dev_combinations, THR, deviceID, deviceCount);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime,start,stop);
	fprintf(stderr,"kernel computation terminated on GPU %d. GPU-Time required (ms): %4.5f\n", deviceID, elapsedTime);
	fprintf(stderr,"GPU %d, elapsed %f seconds\n", deviceID, ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
	HANDLE_ERROR( cudaEventDestroy( start ) );
	HANDLE_ERROR( cudaEventDestroy( stop ) );
	
	
	cudaMemcpy(output + s, dev_output, d * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(tp + s, dev_tp, d * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(fp + s, dev_fp, d * sizeof(int), cudaMemcpyDeviceToHost);
	
	//fprintf(stderr,"GPU%d, output, tp, fp copied to host \n", deviceID);
	
	
  	//free
  	cudaFree(dev_mat_SNP);
	//cudaFree(dev_v_pheno);
	cudaFree(dev_output);
	cudaFree(dev_tp);
	cudaFree(dev_fp);
	cudaFree(dev_combinations);
	//cudaFree(dev_cv_indices);
	//printf("GPU %d, returning; \nelapsed %f seconds\n", deviceID, ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
	return 0;
}

/************************/
//MAIN
/************************/

int main(int argc, char **argv)
{

	clock_t start_clock = clock();
	
	
	//Parsing the input parameters
	parseArgs(argc,argv);
	 
	 
  	//print_cudaGetDeviceProperties(); 
  	
  	int deviceCount;
	HANDLE_ERROR( cudaGetDeviceCount( &deviceCount ) );
	//printf( "found %d devices\n", deviceCount );
	if (NUMDEVICES > deviceCount || deviceCount == 0){
		fprintf(stderr,"error: less devices detected (%d) than specified (%d)! Exiting...\n", deviceCount, NUMDEVICES);
		return 0;
	}
	//fprintf(stderr,"Gonna use first %d of %d devices\n", NUMDEVICES, deviceCount);
	deviceCount = imin(deviceCount, NUMDEVICES);
  	
  	cudaDeviceProp prop;
  	HANDLE_ERROR( cudaGetDeviceProperties( &prop, 0 ) );
  	
  	if (NUMCOMBS/deviceCount > prop.maxGridSize[0])
  		fprintf(stderr,"Supported up to %d*%d combs. input'll be considered up to that combination. Run again with new file later.\n", deviceCount, prop.maxGridSize[0]);
  	
  	
  	//fprintf(stderr,"\n*****************\n");
	//fprintf(stderr,"Multifactor Dimensionality Reduction\n");
	//fprintf(stderr,"*****************\n\n");	
	
	//Allocate host memory 
	int* mat_SNP = (int*)malloc(mat_SNP_size); 
	int* v_pheno = (int*)malloc(v_pheno_size);
	float* output = (float*)malloc(output_size);
	int* tp = (int*)malloc(tp_size);
	int* fp = (int*)malloc(fp_size);
	int* combinations = (int*)malloc(combinations_size);
	int* cv_indices = (int*)malloc(indices_size);
	
	/*
	if (CV <= 1){
		fprintf(stderr,"will run only one pass, no train-test... \n");
	}
	*/
	//generate a permutation of the individuals indices
	for(int i=0;i<NIND;++i){
        	*(cv_indices + i) = i;
    		}
    	
    	/*	
		//permute r with Fisher-Yates shuffling algorithm
	for (int i = NIND; i >= 0; --i){
		//generate a random number [0, n-1]
		int j = rand() % (i+1);

		//swap the last element with element at random index
		int temp = *(cv_indices + i);
		*(cv_indices + i) = *(cv_indices + j);
		*(cv_indices + j) = temp;
	}
	*/
  	
  	//Read the matrix in host data
	readintData(genoFile, NSNPS, NIND, mat_SNP);
	//fprintf(stderr,"geno file read..\n");
	readintData(phenoFile, NIND, 1, v_pheno);
	//fprintf(stderr,"pheno file read..\n");
	
	int ncases;
	if (THR < 0){
		ncases = 0;
		for (int i = 0; i < NIND; i++){
			if ( *(v_pheno + i) )
				ncases += 1;
		}
		THR = float(ncases)/(NIND - ncases);
		printf("no input threshold; automatically set to %f = %d/%d \n", THR, ncases, NIND-ncases);
	}
	
	
	if (basic_model % 2 != 0){
		
		//fprintf(stderr,"detected NULL model. Exhaustive in-memory generation...\n");
		int t = 0;
		if (ORDER == 2) {
			for(int i=0;i<NSNPS;++i){
				for(int j=i+1;j<NSNPS;++j){
					combinations[t] = i;
					combinations[t+1] = j;
					t += ORDER;
					
				}
			}
		}
		else if (ORDER == 3) {
			for(int i=0;i<NSNPS;++i){
				for(int j=i+1;j<NSNPS;++j){
					for(int k=j+1;k<NSNPS;++k){
							combinations[t] = i;
							combinations[t+1] = j;
							combinations[t+2] = k;
							t += ORDER;
					}
				}
			}
		
		}
		else if (ORDER == 4) {
			for(int i=0;i<NSNPS;++i){
				for(int j=i+1;j<NSNPS;++j){
					for(int k=j+1;k<NSNPS;++k){
						for(int l=k+1;l<NSNPS;++l){
							combinations[t] = i;
							combinations[t+1] = j;
							combinations[t+2] = k;
							combinations[t+3] = l;
							t += ORDER;
						}
					}
				}
			}
		}
		else if (ORDER == 5) {
			for(int i=0;i<NSNPS;++i){
				for(int j=i+1;j<NSNPS;++j){
					for(int k=j+1;k<NSNPS;++k){
						for(int l=k+1;l<NSNPS;++l){
							for(int m=l+1;m<NSNPS;++m){
								combinations[t] = i;
								combinations[t+1] = j;
								combinations[t+2] = k;
								combinations[t+3] = l;
								combinations[t+4] = m;
								t += ORDER;
							}
						}
					}
				}
			}
		}
		else if (ORDER == 6) {
			for(int i=0;i<NSNPS;++i){
				for(int j=i+1;j<NSNPS;++j){
					for(int k=j+1;k<NSNPS;++k){
						for(int l=k+1;l<NSNPS;++l){
							for(int m=l+1;m<NSNPS;++m){
								for(int n=m+1;n<NSNPS;++n){
									combinations[t] = i;
									combinations[t+1] = j;
									combinations[t+2] = k;
									combinations[t+3] = l;
									combinations[t+4] = m;
									combinations[t+5] = m;
									t += ORDER;
								}
							}
						}
					}
				}
			}
		}
	//fprintf(stderr,"combinations in-memory generated..\n");	
	}
	else{	
		
		//Read combinations
		fprintf(stderr,"detected disease model\n");
		readCombinations(combFile, NUMCOMBS, ORDER, combinations);
		fprintf(stderr,"combinations file read..\n");
	}
	

			
  	//printf("calling therads; \nelapsed %f seconds\n", ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
  	CUTThread threads[deviceCount-1];
  	DataStruct  data[deviceCount];
  	
  	for(int i=0;i<deviceCount;++i){
  		data[i].deviceID = i;
  		data[i].deviceCount = deviceCount;
		data[i].mat_SNP = mat_SNP;
		data[i].combinations = combinations; //same for all. split happens in the CUDAMemCpy
		data[i].v_pheno = v_pheno;
		data[i].cv_indices = cv_indices;
		data[i].output = output;
		data[i].tp = tp;
		data[i].fp = fp;
		data[i].start_clock = start_clock;
  		if (i != (deviceCount -1)){
  			threads[i] = start_thread( routine, &(data[i]) );
  			
  		}
  		else{
  			routine( &(data[i]) );
  			for(int j=0;j<deviceCount-1;++j)
				end_thread( threads[j] );
		}
  	}
  	
	//printf("last thread terminated. \nelapsed %f seconds\n", ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
	free(mat_SNP);
	free(v_pheno);
	free(cv_indices);
	
	
	if (strcmp(outputFile, "no_out")){
		//sort output and print to file
		struct str **array = (struct str **) malloc(NUMCOMBS * sizeof(struct str*));
		
	
		FILE *fpout;
		fpout = fopen(outputFile, "w");
	
		if (CV == 1){
			//fprintf(fpout,"---------- measure ----------\n");
			
			for(int i=0;i<NUMCOMBS;i++){
				*(array + i) = (struct str*) malloc(sizeof(struct str));
				(*(array + i))->value = *(output + i);
				(*(array + i))->index = i;
				//objects[i].value=*(output + i);
				//objects[i].index=i;
			}
			//printf("last thread terminated. \nelapsed %f seconds\n", ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
			fprintf(stderr,"Effectively sorting output..elapsed %f seconds\n", ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
			

			//fprintf(stderr,"\nelapsed %f seconds\n", ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
		       	//fprintf(stderr,"Effectively writing to file: %s \n", outputFile);
		       	int cut = imin(CUT,NUMCOMBS);

			/*
			//version merge all

		       	mergeSort(array, 0, NUMCOMBS - 1);	
		       	fprintf(stderr,"terminated merge sort.. elapsed %f seconds\n", ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
		       	//fprintf(stderr,"breaking to first %d combinations", cut);
			for (int j = 0; j < cut; j++){
				for (int q=0; q< ORDER; q++){
					if (q == 0)
						fprintf(fpout,"snp%d ", *(combinations + (*(array + j))->index * ORDER + q));
					else if (q == ORDER -1){
					
						int tp_val = *(tp + (*(array + j))->index);
						int fp_val = *(fp + (*(array + j))->index);
						int tn_val = NIND - ncases - fp_val;
						int fn_val = ncases - tp_val;
						
						fprintf(fpout,"snp%d %f %d %d %d %d\n", 
							//*(combinations + j * ORDER + q), *(output + NUMCOMBS * (CV -1) + 1 * objects[j].index + 0));
							//*(combinations + j * ORDER + q), objects[j].value);
							*(combinations +(*(array + j))->index * ORDER + q), (*(array + j))->value, tp_val, fp_val, tn_val, fn_val);
					}
					else
						fprintf(fpout,"snp%d ", *(combinations + (*(array + j))->index * ORDER + q));
					
			
				}
			}
			//end version merge all
			
			*/
			
			
			
			//alternative
			//new version
		       	struct str **top_cut = (struct str **) malloc(cut * sizeof(struct str*));
		       	for (int j = 0; j < cut; j++){
		       		*(top_cut + j) = (struct str*) malloc(sizeof(struct str));
		       		(*(top_cut + j))->value = 0;
		       	}
		       	
		       	//scan through output and insert in top_cut
		       	float min_cut = 0;
		       	int ind;
		       	float* a = (float*) malloc(2 * sizeof(float));
		       	
		       	for (int j = 0; j < NUMCOMBS; j++){
		       		
		       		
		       		if (*(output + j) < min_cut)
		       			continue;
		       		//fprintf(stderr,"---\n");
		       		//fprintf(stderr,"min_cut,index: %f %d\n", min_cut, ind);
		       		else {
		       			
		       			extract_min(top_cut, cut, a);
			       		//fprintf(stderr,"%f, %f", *(a+0), *(a+1));
			       		min_cut =*(a+0);
			       		ind = (int) *(a+1);
		       			
		       			//fprintf(stderr,"removing %f %d\n,", min_cut, ind);
		       			(*(top_cut + ind))-> index = j;
		       			(*(top_cut + ind))-> value = *(output + j);
		       			//fprintf(stderr,"inserting %f %d",*(output + j), j);
		       			//for (int r = 0; r < cut; r++){
		       				//fprintf(stderr,"index r-th,value r-th: %d,%f \n", (*(top_cut + r))-> index, (*(top_cut + r))-> value);
		       			//fprintf(stderr,"\n");
		       			//}
		       			//fprintf(stderr,"inserted output + j at index j: %f, %d\n", *(output + j), j);
		       		}
		       	}

		       	fprintf(stderr,"top cut sorting... \n");
		       	
		       mergeSort(top_cut, 0, cut - 1);	
			
			for (int j = 0; j < cut; j++)
		       		fprintf(stderr,"index,value: %d,%f \n", (*(top_cut + j))-> index, (*(top_cut + j))-> value);
		       fprintf(stderr,"top cut sorted; elapsed %f seconds\n", ((float)(clock() - start_clock) / CLOCKS_PER_SEC));
		       
		       
		       
		       for (int j = 0; j < cut; j++){
				for (int q=0; q< ORDER; q++){
					if (q == 0)
						fprintf(fpout,"snp%d ", *(combinations + (*(top_cut + j))->index * ORDER + q));
					else if (q == ORDER -1){
					
						int tp_val = *(tp + (*(top_cut + j))->index);
						int fp_val = *(fp + (*(top_cut + j))->index);
						int tn_val = NIND - ncases - fp_val;
						int fn_val = ncases - tp_val;
						
						fprintf(fpout,"snp%d %f %d %d %d %d\n", 
							//*(combinations + j * ORDER + q), *(output + NUMCOMBS * (CV -1) + 1 * objects[j].index + 0));
							//*(combinations + j * ORDER + q), objects[j].value);
							*(combinations +(*(top_cut + j))->index * ORDER + q), (*(top_cut + j))->value, tp_val, fp_val, tn_val, fn_val);
					}
					else
						fprintf(fpout,"snp%d ", *(combinations + (*(top_cut + j))->index * ORDER + q));
					
			
				}
			}
			
			//end new version
			
			
		
		/*
		else{
	  		for (int cv = 0; cv < CV; cv++){
		  		fprintf(fpout,"---------- CV %d/%d train_measure test_measure(s) ----------\n", cv+1, CV);
		  		
		  		for(int i=0;i<NUMCOMBS;i++){
					objects[i].value=*(output + NUMCOMBS * cv + 2 * i + 1); //sorting on test measure!
					objects[i].index=i;
				}
				qsort(objects,NUMCOMBS,sizeof(objects[0]),cmp);

		  		for (int j = 0; j < NUMCOMBS; j++){
					for (int q=0; q< ORDER; q++){
						if (q == 0)
							fprintf(fpout,"snp%d ", *(combinations + j * ORDER + q));
						else if (q == ORDER -1)
							fprintf(fpout,"snp%d %f %f\n", 
								*(combinations + j * ORDER + q),
								*(output +  NUMCOMBS * cv + 2 * objects[j].index + 0), //sorted on test measure!
								*(output + NUMCOMBS * cv + 2 * objects[j].index + 1) );
								//*(output +  NUMCOMBS * cv + 2 * j + 0),
								//*(output + NUMCOMBS * cv + 2 * j + 1) );
						else
							fprintf(fpout,"snp%d ", *(combinations + j * ORDER + q));
					}

				}
			}
		}*/
		//fprintf(stderr,"Output written to file %s\n", outputFile);
		}
	}
	
	//else
		//fprintf(stderr,"Output was not saved to file \n");
	
	free(output);
	free(tp);
	free(fp);
   	free(combinations);

 	return 0;
}
