/* File that contains the variable declarations */
#include "geras.h" 
#include <stdio.h>

/* First function called. It initiallizes all the functions and variables used by GERAS */
void geras_init(int geras, int start_search){
        int i, startThreads=0;
        int numCores = sysconf(_SC_NPROCESSORS_ONLN);
        /*Initialization of RAPL */
        //geras_detect_cpu(); 
        geras_detect_packages();
        /*End initialization of RAPL */

        /* Initialization of the variables necessary to perform the search algorithm */
		if(start_search == 0){
		        startThreads = numCores;
		        while(startThreads != 2 && startThreads != 3 && startThreads != 5){
		                startThreads = startThreads/2;
		        }
		}else{
			startThreads = start_search;
		}
        for(i=0;i<MAX_KERNEL;i++){
                gerasKernels[i].startThreads = startThreads;
                gerasKernels[i].numThreads = numCores;
                gerasKernels[i].numCores = numCores;
                gerasKernels[i].initResult = 0.0;
                gerasKernels[i].state = REPEAT;
                gerasKernels[i].steps = 0;
                gerasKernels[i].bestTime = 1000.00;
		gerasKernels[i].total_region_perf = 0.0;
                gerasKernels[i].gerasMetric = geras;
                gerasKernels[i].lastResult = 0.0;
                idKernels[i]=0;
        }

        /* Start the counters for energy and time for all the application execution */
        id_actual_region = MAX_KERNEL-1;
        geras_start_amd_msr();
        initGlobalTime = omp_get_wtime();

}

/* It defines the number of threads that will execute the actual parallel region based on the current state of the search algorithm */
int geras_resolve_num_threads(uintptr_t ptr_region){
        int i;
        id_actual_region = -1;
        /* Find the actual parallel region */
        for(i=0;i<totalKernels;i++){
                if(idKernels[i] == ptr_region){
                        id_actual_region=i;
                        break;
                }
        }
        /* If a new parallel region is discovered */
        if(id_actual_region == -1){
                idKernels[totalKernels] = ptr_region;
                id_actual_region = totalKernels;
                totalKernels++;
        }
        /* Check the state of the search algorithm. */
        switch(gerasKernels[id_actual_region].state){
                case END:
                        gerasKernels[id_actual_region].initResult = omp_get_wtime();  /* It is useful only if the continuous adaptation is enable. Otherwise, it can be disabled */
                        return gerasKernels[id_actual_region].bestThread;
                default:
                        geras_start_amd_msr();
                        gerasKernels[id_actual_region].initResult = omp_get_wtime();
                        return gerasKernels[id_actual_region].numThreads;
        }

}

/* It is responsible for performing the search algorithm */
void geras_end_parallel_region(){
        double time, energy, result=0, ratio;
        if(gerasKernels[id_actual_region].state !=END){
                /* Check the metric that is being evaluated and collect the results */
                switch(gerasKernels[id_actual_region].gerasMetric){
                        case PERFORMANCE:
                                time = omp_get_wtime();
                                result = time - gerasKernels[id_actual_region].initResult;
                                break;
                        case AGING:
                        		time = omp_get_wtime() - gerasKernels[id_actual_region].initResult;
                        		energy = geras_end_amd_msr();
					result = sqrt( (energy*energy)+(time*time));
					gerasKernels[id_actual_region].total_region_perf += time;
	                                gerasKernels[id_actual_region].steps++;
        	                        if(time < gerasKernels[id_actual_region].bestTime)
                	                        gerasKernels[id_actual_region].bestTime = time;
                        		if(result == 0.00000 || result < 0){
                                        	gerasKernels[id_actual_region].state = REPEAT;
                                        	gerasKernels[id_actual_region].gerasMetric = PERFORMANCE;
                                	}
                                break;
                }
                switch(gerasKernels[id_actual_region].state){
                        case REPEAT:
                                gerasKernels[id_actual_region].state = S0;
                                gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].startThreads; 
                                gerasKernels[id_actual_region].lastThread = gerasKernels[id_actual_region].numThreads;
                                break;
                        case S0:
                                gerasKernels[id_actual_region].bestResult = result;
                                gerasKernels[id_actual_region].bestThread = gerasKernels[id_actual_region].numThreads;
                                gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].bestThread*2; 
                                gerasKernels[id_actual_region].state = S1;
                                break;
                        case S1:
                                if(result < gerasKernels[id_actual_region].bestResult){ //comparing S0 to REPEAT
                                        gerasKernels[id_actual_region].bestResult = result;
                                        gerasKernels[id_actual_region].bestThread = gerasKernels[id_actual_region].numThreads;
                                        /* if there are opportunities for improvements, then double the number of threads */
                                        if(gerasKernels[id_actual_region].numThreads * 2 <= gerasKernels[id_actual_region].numCores){
                                                gerasKernels[id_actual_region].lastThread = gerasKernels[id_actual_region].numThreads;
                                                gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].bestThread*2;
                                                gerasKernels[id_actual_region].state = S1;
                                        }else{
                                                /* It means that the best number so far is equal to the number of cores */
                                                /* Then, it will realize a guided search near to this number */
                                                gerasKernels[id_actual_region].pass = gerasKernels[id_actual_region].lastThread/2;                                                
                                                gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].numThreads  - gerasKernels[id_actual_region].pass;
                                                if(gerasKernels[id_actual_region].pass == 1)
                                                        gerasKernels[id_actual_region].state = S3;        
						else
	                                                gerasKernels[id_actual_region].state = S2;
                                        }
                                }else{ 
                                        /* Thread scalability stopped */
                                        /* Find the interval of threads that provided this result. */
                                        /* if the best number of threads so far is equal to the number of cores, then go to.. */
                                        if(gerasKernels[id_actual_region].bestThread == gerasKernels[id_actual_region].numCores/2){
                                                gerasKernels[id_actual_region].pass = gerasKernels[id_actual_region].lastThread/2;                                                
                                                gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].numThreads  - gerasKernels[id_actual_region].pass;
                                                if(gerasKernels[id_actual_region].pass == 1)
                                                        gerasKernels[id_actual_region].state = S3;        
        					else
		                                        gerasKernels[id_actual_region].state = S2;
                                        }else{
                                                gerasKernels[id_actual_region].pass = gerasKernels[id_actual_region].lastThread/2;                                                
                                                gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].numThreads  + gerasKernels[id_actual_region].pass;
                                                if(gerasKernels[id_actual_region].pass == 1)
                                                        gerasKernels[id_actual_region].state = S3;        
                                                else
														gerasKernels[id_actual_region].state = S2;  
                                        }

                                }
                                break;
                        case S2:        
                                if(gerasKernels[id_actual_region].bestResult < result){
                                        gerasKernels[id_actual_region].pass = gerasKernels[id_actual_region].pass/2;
                                        gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].numThreads + gerasKernels[id_actual_region].pass;
                                        if(gerasKernels[id_actual_region].pass == 1)
                                                gerasKernels[id_actual_region].state = S3;
                                        else
                                                gerasKernels[id_actual_region].state = S2;
                                }else{
                                        gerasKernels[id_actual_region].bestThread = gerasKernels[id_actual_region].numThreads;
                                        gerasKernels[id_actual_region].bestResult = result;
                                        gerasKernels[id_actual_region].pass = gerasKernels[id_actual_region].pass/2;
                                        gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].numThreads + gerasKernels[id_actual_region].pass;
                                        if(gerasKernels[id_actual_region].pass == 1)
                                                gerasKernels[id_actual_region].state = S3;
                                        else
                                                gerasKernels[id_actual_region].state = S2;
                                }
                                break;                        
                        case S3: //The last comparison to define the best number of threads
                                if(result < gerasKernels[id_actual_region].bestResult){
                                        gerasKernels[id_actual_region].bestThread = gerasKernels[id_actual_region].numThreads;
                                        gerasKernels[id_actual_region].state = END;
                                }else{
                                        gerasKernels[id_actual_region].state = END;
                                }
                                break;
                }
        }else{
                result = omp_get_wtime() - gerasKernels[id_actual_region].initResult;
                if(result > 0.1){
                        if(gerasKernels[id_actual_region].lastResult == 0.0){
                                gerasKernels[id_actual_region].lastResult = result;
                        }else{
                                ratio = gerasKernels[id_actual_region].lastResult/result;
                                if(ratio < 0.7 || ratio > 1.3){
                                        gerasKernels[id_actual_region].numThreads = gerasKernels[id_actual_region].numCores;
                                        gerasKernels[id_actual_region].lastResult = 0.0;
                                        gerasKernels[id_actual_region].state = REPEAT;
                                }else{
                                        gerasKernels[id_actual_region].lastResult = result;
                                }
                        }
                }
        }
}

/* It finalizes the environment of Aurora */
void geras_destructor(){
        double time = omp_get_wtime() - initGlobalTime;
        id_actual_region = MAX_KERNEL-1;
        double energy = geras_end_amd_msr();
        printf("GERAS - Execution Time: %.5f seconds\n", time);
        printf("GERAS - Energy: %.5f joules\n",energy);
}


void geras_detect_packages() {

	char filename[STRING_BUFFER];
	FILE *fff;
	int package;
	int i;

	for(i=0;i<MAX_PACKAGES;i++) package_map[i]=-1;

	for(i=0;i<MAX_CPUS;i++) {
		sprintf(filename,"/sys/devices/system/cpu/cpu%d/topology/physical_package_id",i);
		fff=fopen(filename,"r");
		if (fff==NULL) break;
		fscanf(fff,"%d",&package);
		fclose(fff);

		if (package_map[package]==-1) {
			gerasTotalPackages++;
			package_map[package]=i;
		}

	}

}

void geras_start_amd_msr(){
	char msr_filename[STRING_BUFFER];
	int fd;
	sprintf(msr_filename, "/dev/cpu/0/msr");
	fd = open(msr_filename, O_RDONLY);
	if ( fd < 0 ) {
		if ( errno == ENXIO ) {
			fprintf(stderr, "rdmsr: No CPU 0\n");
			exit(2);
		} else if ( errno == EIO ) {
			fprintf(stderr, "rdmsr: CPU 0 doesn't support MSRs\n");
			exit(3);
		} else {
			perror("rdmsr:open");
			fprintf(stderr,"Trying to open %s\n",msr_filename);
			exit(127);
		}
	}
	uint64_t data;
	pread(fd, &data, sizeof data, AMD_MSR_PACKAGE_ENERGY);
	//gerasKernels[id_actual_region].kernelBefore[0] = read_msr(fd, AMD_MSR_PACKAGE_ENERGY);
	gerasKernels[id_actual_region].kernelBefore[0] = (long long) data;
}

double geras_end_amd_msr(){
	char msr_filename[STRING_BUFFER];
        int fd;
        sprintf(msr_filename, "/dev/cpu/0/msr");
        fd = open(msr_filename, O_RDONLY);
		uint64_t data;
        pread(fd, &data, sizeof data, AMD_MSR_PWR_UNIT);
	int core_energy_units = (long long) data;
	unsigned int energy_unit = (core_energy_units & AMD_ENERGY_UNIT_MASK) >> 8;
	pread(fd, &data, sizeof data, AMD_MSR_PACKAGE_ENERGY);
	gerasKernels[id_actual_region].kernelAfter[0] = (long long) data;
	double result = (gerasKernels[id_actual_region].kernelAfter[0] - gerasKernels[id_actual_region].kernelBefore[0])*pow(0.5,(float)(energy_unit));
	return result;
}
