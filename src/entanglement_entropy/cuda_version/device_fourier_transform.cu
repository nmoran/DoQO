#include "device_fourier_transform.h"
#include "host_fourier_transform.h"
 
#define CUDA_BLOCKS 510  		//this should be a multiple of 30 
#define CUDA_THREADS_PER_BLOCK 256	//this should be a multiple of 32 
#define MAX_FILLING 10
#define MAX_SITES 20

#define PI 3.141592653589793f

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
__device__ void get_monomial_bosons_cuda(unsigned long Index, int Sites, int Bosons, int *Monomial) 
  {
    int mom = 1;
    int idx = 0;
    int pos = 0;
    while ( pos < (Sites + Bosons) ) {
      if ( (Index & (1ul << pos) ) > 0ul ) {
	  Monomial[idx++] = mom;
      } else {
	  mom += 1;
      }	  
      pos++;
    }  
  }

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order and also return the normal.
__device__ void get_monomial_bosons_cuda(unsigned long Index, int Sites, int Bosons, int *Monomial, float *Normal) 
  {
    int mom = 1;
    int idx = 0;
    int pos = 0;
    *Normal = 1.0f;
    float num = 0;
    while ( pos < (Sites + Bosons) ) {
      if ( (Index & (1ul << pos) ) > 0ul ) {
	  Monomial[idx++] = mom;
	  num += 1.0f;
	  *Normal *= num;
      } else {
	  mom += 1;
	  //if ( num > 1.0 ) *normal *= num; 
	  num = 0.0f; 
      }	  
      pos++;
    }
    //if ( num > 1.0 ) *normal *= num; 
    *Normal = sqrtf(*Normal);
  }

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order and also return the normal.
unsigned long convert_from_monomial_bosons(int Sites, int Bosons, int *Monomial) 
  {
    unsigned long config = 0;
    
    for ( int i = 0 ; i < Bosons ; i++ ) 
      {
	config += 1ul << (i + Monomial[Bosons-i-1]);
      }
    return config;     
  }


// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order and also return the normal.
__device__ unsigned long convert_from_monomial_bosons_cuda(int Sites, int Bosons, int *Monomial) 
  {
    unsigned long config = 0;
    
    for ( int i = 0 ; i < Bosons ; i++ ) 
      {
	config += 1ul << (i + Monomial[Bosons-i-1]);
      }
    return config;     
  }


// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
__device__ void get_monomial_fermions_cuda(unsigned long Index, int Sites, int Fermions, int *Monomial) 
  {
    int mom = 1;
    int idx = 0;
    int pos = 0;
    while ( pos < Sites ) {
      if ( (Index & (1ul << pos) ) > 0ul ) {
	  Monomial[idx++] = mom;
	  mom += 1; 
      } else {
	  mom += 1;
      }	  
      pos++;
    }  
  }

/*__device__ int permute_cuda(int *Monomial, int Len)
  {
    int k,l,tmp;
    
    k = Len - 2;
    while ((k > 0) && (Monomial[k] >= Monomial[k+1]) ) k--;
    
    if ( (k == 0) && (Monomial[0] >= Monomial[1]) ) return 0;
    
    l = Len - 1;
    while ( Monomial[l] <= Monomial[k] ) l--;
    
    tmp = Monomial[k];
    Monomial[k] = Monomial[l];
    Monomial[l] = tmp;
    
    l = Len - 1;
    k = k + 1;
    
    while ( k < l ) 
    {
      tmp = Monomial[k];
      Monomial[k] = Monomial[l];
      Monomial[l] = tmp;
      l--;
      k++;
    }    
    return 1;       
  }*/
  
__device__ int permute_cuda(int *Monomial, int Len)
  {
    int k,l,tmp;
    k = Len - 2;
    while ((k > 0) && (Monomial[k] >= Monomial[k+1]) ) k--;
    if ( (k == 0) && (Monomial[0] >= Monomial[1]) ) return 0;
    l = Len - 1;
    while ( Monomial[l] <= Monomial[k] ) l--;
    tmp = Monomial[k];
    Monomial[k] = Monomial[l];
    Monomial[l] = tmp;
    l = Len - 1;
    k = k + 1;
    while ( k < l ) 
    {
      tmp = Monomial[k];
      Monomial[k] = Monomial[l];
      Monomial[l] = tmp;
      l--;
      k++;
    }    
    return 1;       
  }  
  
__device__ float configuration_normal_cuda(unsigned long Config, int NumberSites, unsigned long Mask)
  {
    float normal = 1.0f;        
    for ( int i = 1 ; i < NumberSites ; i++ ) 
    {
      if ( Config == ( ((Config << i) | (Config >> (NumberSites - i))) & Mask ) ) 
	{
	  normal += 1.0f;
	}
    }
    return (1.0f/sqrtf(normal));
  }


__global__ void test_monomial_conversion(int *fmonomial, int NumberSites, int Filling, unsigned long Config)
  {
//    float stateNormal;
    get_monomial_fermions_cuda(Config, NumberSites, Filling, fmonomial);
    //while ( permute_cuda(fmonomial, Filling) ) {}     
  }

__global__ void calculate_FT_elements_alt(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, float *PhaseArrayR, int offset, float pre_norm) 
  {	    
    /*float l_PhaseArrayR[MAX_SITES];
    for ( int i = 0 ; i < NumberSites ; i++ ) 
      {
	l_PhaseArrayR[i] = PhaseArrayR[i];
      }*/
    int fMonomialPtr[MAX_FILLING];
    int fMonomialTmpPtr[MAX_FILLING];
    int oMonomialPtr[MAX_FILLING];
    __shared__ float shared_tmp_vals[CUDA_THREADS_PER_BLOCK];
    int fidx = blockIdx.x + offset * gridDim.x;
    if ( fidx < FTBasisDim ) 
      {	
	shared_tmp_vals[threadIdx.x] = 0.0f;
	float stateNormal;
	get_monomial_bosons_cuda(FTBasisElements[fidx], NumberSites, Filling, fMonomialPtr);
	get_monomial_bosons_cuda(FTBasisElements[fidx], NumberSites, Filling, fMonomialTmpPtr, &stateNormal);
	
	int oidx = threadIdx.x;
	while (oidx < VectorLength )
	  {
	    int phaseIdx = 0;
	    get_monomial_fermions_cuda(BasisElements[oidx], NumberSites, Filling, oMonomialPtr);	    
            for ( int i = 0; i < Filling; i++ )
              {
	        fMonomialPtr[i] = fMonomialTmpPtr[i];
		phaseIdx += oMonomialPtr[i] * fMonomialPtr[i];
              }
            float tmp_val2 = 0.0f;
	    int permuteFlg = 1, k, l, tmp;
	    while ( permuteFlg ) 
	      {
		tmp_val2 += PhaseArrayR[phaseIdx % NumberSites];
    
		k = Filling - 2;
		while ((k > 0) && (fMonomialPtr[k] >= fMonomialPtr[k+1]) ) k--;
		
		if ( (k == 0) && (fMonomialPtr[0] >= fMonomialPtr[1]) ) 
		  {
		    permuteFlg = 0;
		  } 
		else 
		  {
		    l = Filling - 1;
		    while ( fMonomialPtr[l] <= fMonomialPtr[k] ) l--;
		  
		    phaseIdx += (fMonomialPtr[l] - fMonomialPtr[k]) * oMonomialPtr[k] + (fMonomialPtr[k] - fMonomialPtr[l]) * oMonomialPtr[l];
		    tmp = fMonomialPtr[k];
		    fMonomialPtr[k] = fMonomialPtr[l];
		    fMonomialPtr[l] = tmp;
		    l = Filling - 1;
		    k = k + 1;
		    
		    while ( k < l ) 
		    {
		      phaseIdx += (fMonomialPtr[l] - fMonomialPtr[k]) * oMonomialPtr[k] + (fMonomialPtr[k] - fMonomialPtr[l]) * oMonomialPtr[l];
		      tmp = fMonomialPtr[k];
		      fMonomialPtr[k] = fMonomialPtr[l];
		      fMonomialPtr[l] = tmp;
		      l--;
		      k++;
		    }    
		    permuteFlg = 1;    
		  }
	      }
	    
	    shared_tmp_vals[threadIdx.x] += tmp_val2 * VectorElements[oidx] ;
	    oidx += blockDim.x;
	  }
    	  
	__syncthreads();
	int i = blockDim.x/2;
	while ( i != 0 ) 
	  {
	    if ( threadIdx.x < i )
	      {
	        shared_tmp_vals[threadIdx.x] += shared_tmp_vals[threadIdx.x + i];
	      }
	    __syncthreads();  
	    i /= 2;
	  }
	if ( threadIdx.x == 0 )  
	  {
	    FTVectorElements[fidx] = shared_tmp_vals[0] * pre_norm * stateNormal;
	  }
      }      
  }
  
__global__ void calculate_FT_elements(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, float *PhaseArrayR, int offset, float pre_norm) 
  {	            
    int fidx;     
    int fMonomialPtr[MAX_FILLING];
    int fMonomialTmpPtr[MAX_FILLING];
    int oMonomialPtr[MAX_FILLING];
    __shared__ float shared_tmp_vals[CUDA_THREADS_PER_BLOCK];
    float l_PhaseArrayR[MAX_SITES];
    
    for ( int i = 0 ; i < NumberSites ; i++ )
      {
        l_PhaseArrayR[i] = PhaseArrayR[i];
      }
    
    int  phaseIdx;    
    float stateNormal, tmp_val2 ;
        
    fidx = blockIdx.x + offset * gridDim.x;
    if ( fidx < FTBasisDim ) 
      {	
	shared_tmp_vals[threadIdx.x] = 0.0f;
	get_monomial_bosons_cuda(FTBasisElements[fidx], NumberSites, Filling, fMonomialPtr, &stateNormal);
	get_monomial_bosons_cuda(FTBasisElements[fidx], NumberSites, Filling, fMonomialTmpPtr, &stateNormal);
	
	int oidx = threadIdx.x;
	while (oidx < VectorLength )
	  {
	    get_monomial_fermions_cuda(BasisElements[oidx], NumberSites, Filling, oMonomialPtr);	    
            for ( int i = 0; i < Filling; i++ )
              {
	        fMonomialPtr[i] = fMonomialTmpPtr[i];
              }
            tmp_val2 = 0.0f;        
            do 
              {
	        phaseIdx = 0;
	        for ( int i = 0 ; i < Filling; i++ )
	          {
	            phaseIdx += oMonomialPtr[i] * fMonomialPtr[i];	            
	          }
	        tmp_val2 += l_PhaseArrayR[phaseIdx % NumberSites] ;
	      }
	    while ( permute_cuda(fMonomialPtr, Filling) );
	    shared_tmp_vals[threadIdx.x] += tmp_val2 * VectorElements[oidx] ;
	    oidx += blockDim.x;
	  }
    	  
	__syncthreads();
	int i = blockDim.x/2;
	while ( i != 0 ) 
	  {
	    if ( threadIdx.x < i )
	      {
	        shared_tmp_vals[threadIdx.x] += shared_tmp_vals[threadIdx.x + i];
	      }
	    __syncthreads();  
	    i /= 2;
	  }
	if ( threadIdx.x == 0 )  
	  {
	    FTVectorElements[fidx] = shared_tmp_vals[0] * pre_norm * stateNormal;
	  }
      }      
  }  
 
void fourier_transform_vector_cuda(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, int CurrentMomentum, int SpinMomentum, bool single, int CudaDevice, string OutputPrefix, int Offset)
  {        
    cudaError_t cudaerr;
    //get device code
//     cudaDeviceProp prop;
//     int dev;
//     cudaGetDevice( &dev );    
//     memset( &prop, 0, sizeof( cudaDeviceProp ) );
//     prop.multiProcessorCount = 16;    
//     cudaChooseDevice( &dev, &prop );
//     cout << "ID of CUDA device with closest to 32 cores. " << dev << endl;    
    cudaerr = cudaSetDevice( CudaDevice );
    cout << "Setting CUDA device to: " << CudaDevice << ": " << cudaGetErrorString(cudaerr) << endl;  
    
    //prepare phase arrays
    float *phaseArrayR = new float[NumberSites+1];
    float *phaseArrayI = new float[NumberSites+1];
    phaseArrayR[0] = 1.0f;
    phaseArrayI[0] = 0.0f;
    for( int i = 1  ; i <= NumberSites ; i++ ) 
      {	    	    				
        float angle = ( PI * 2.0f * (float)i) / ((float)NumberSites);
	phaseArrayR[i] = (float)cos(angle);
	phaseArrayI[i] = (float)sin(angle);	
      }
      
    float norm, momentumPhase = 1.0f;    
    if ( SpinMomentum != -1 ) 
      {
	for ( int i = 1; i < NumberSites ; i++ )
	  {
	    momentumPhase += phaseArrayR[(i * (CurrentMomentum + SpinMomentum)) % NumberSites];
	  }
      }
    norm = sqrtf(momentumPhase) * sqrtf(1.0f/(float)powf((float)NumberSites,(float)Filling));  
                    
    //declare and allocate space.
    float  *d_FTVectorElements, *d_VectorElements, *d_phaseArrayR;
    unsigned long *d_FTBasisElements, *d_BasisElements;    
            
    cudaerr = cudaMalloc( (void**) &d_phaseArrayR, (NumberSites + 1) * sizeof(float)); cout << "Allocating d_phaseArrayR: " << cudaGetErrorString(cudaerr) << endl;   
    cudaerr = cudaMalloc( (void**) &d_FTBasisElements, FTBasisDim * sizeof(unsigned long) ); cout << "Allocating d_FTBasisELements: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMalloc( (void**) &d_FTVectorElements, FTBasisDim * sizeof(float) );cout << "Allocating d_FTVectorElements: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMalloc( (void**) &d_BasisElements, VectorLength * sizeof(unsigned long) );cout << "Allocating d_BasisElements: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMalloc( (void**) &d_VectorElements, VectorLength * sizeof(float) );cout << "Allocating d_VectorElements: " << cudaGetErrorString(cudaerr) << endl;
    
    cudaerr = cudaMemcpy( d_phaseArrayR, phaseArrayR, (NumberSites + 1) * sizeof(float), cudaMemcpyHostToDevice ); cout << "Copying phaseArrayR: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMemcpy( d_FTBasisElements, FTBasisElements, FTBasisDim * sizeof(unsigned long), cudaMemcpyHostToDevice ); cout << "Copying FTBasisElements: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMemcpy( d_FTVectorElements, FTVectorElements, FTBasisDim * sizeof(float), cudaMemcpyHostToDevice ); cout << "Copying FTVectorElements: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMemcpy( d_BasisElements, BasisElements, VectorLength * sizeof(unsigned long), cudaMemcpyHostToDevice ); cout << "Copying BasisELements: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMemcpy( d_VectorElements, VectorElements, VectorLength * sizeof(float), cudaMemcpyHostToDevice ); cout << "Copying VectorElements: " << cudaGetErrorString(cudaerr) << endl;
	
    
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    int maxOffset = ceil((float)FTBasisDim / (float)CUDA_BLOCKS);
    if ( single ) maxOffset = 1;
    for ( int offset = Offset ; offset < maxOffset ; offset++ ) 
      {
	cout << "Offset: " << offset << " out of " << maxOffset << endl; 
	cudaEventRecord( start, 0 );
	calculate_FT_elements_alt<<<CUDA_BLOCKS, CUDA_THREADS_PER_BLOCK>>>(VectorLength, d_VectorElements, d_BasisElements, FTBasisDim, d_FTVectorElements, d_FTBasisElements, NumberSites, Filling, d_phaseArrayR, offset, norm);
	cudaEventRecord( stop, 0 );
	cudaerr = cudaGetLastError(); cout << "Error after kernel launch: " << cudaGetErrorString(cudaerr) << endl;       
	cudaThreadSynchronize();
	cudaerr = cudaGetLastError(); cout << "Error after kernel finish: " << cudaGetErrorString(cudaerr) << endl;       
	cudaEventSynchronize( stop );
	float elapsedTime;
	cudaEventElapsedTime( &elapsedTime, start, stop);
	cout << "Time from CUDA timer: " << elapsedTime << endl;	
	cudaerr = cudaMemcpy( FTVectorElements, d_FTVectorElements, FTBasisDim * sizeof(float), cudaMemcpyDeviceToHost ); cout << "Copying back FTVectorELements: " << cudaGetErrorString(cudaerr) << endl;       
	string outputFilename;
	stringstream ss("");
	ss << OutputPrefix << "_FT_M_" << CurrentMomentum << "_offset_" << offset << ".vec" ; 
	ofstream out(ss.str().c_str());
	if ( out.is_open() )
	  {
	    out.precision(14);
	    for ( int i = FTBasisDim - 1 ; i >= 0 ; i-- ) 
	      {
		out << scientific << FTVectorElements[i] << endl;
	      }
	    out.close();
	  }	  
      }        
    

    /*for ( int i = 0 ; i < VectorLength ; i++ ) 
      {	
	test_monomial_conversion<<<1,1>>>(d_omonomial, NumberSites, Filling, BasisElements[i]);	    	
	err = cudaMemcpy( omonomial, d_omonomial, Filling * sizeof(int), cudaMemcpyDeviceToHost ); cout << "Copying back fmonomial: " << err << endl;
	
	cout << BasisElements[i] << ": " ;
	for ( int i = 0 ; i < Filling ; i++ ) 
	  {
	    cout << omonomial[i] << " "  ;
	  }
	cout << endl;
      }*/
    cudaEventDestroy( start );
    cudaEventDestroy( stop );
      
        
    cudaFree(d_phaseArrayR);
    cudaFree(d_FTVectorElements);
    cudaFree(d_VectorElements);
    cudaFree(d_FTBasisElements);
    cudaFree(d_BasisElements);
          
    delete [] phaseArrayR;
    delete [] phaseArrayI;
  }

/*__device__ int compare(int *MonomialA, int *MonomialB, int Filling)
  {
    for ( int i = 0 ; i < Filling ; i++ ) 
      {
	if (MonomialA[i] < MonomialB[i] ) 
	  {
	    return 0;
	  } 
	else if (MonomialA[i] > MonomialB[i] ) 
	  {
	    return 1;
	  }
      }
    return 1;
  }
  
  
__global__ void calculate_FT_element(float *VectorElements, unsigned long *BasisElements, int VectorLength, unsigned long FTBasisElement, int NumberSites, int Filling, int CurrentMomentum, int SpinMomentum, float *PhaseArrayR, float *BlockSums, int *Permutations) 
  {	                    
    if ( blockIdx.x < VectorLength ) 
      {
	int fMonomialCurrPtr[MAX_FILLING];
	int fMonomialNextPtr[MAX_FILLING];
	int oMonomialPtr[MAX_FILLING];
	__shared__ float shared_tmp_vals[CUDA_THREADS_PER_BLOCK];
	shared_tmp_vals[threadIdx.x] = 0.0f;    
	
	int  phaseIdx;    
	float stateNormal, norm, momentumPhase;
	unsigned long mask = (1ul << NumberSites) - 1ul;            
	    
	get_monomial_bosons_cuda(FTBasisElement, NumberSites, Filling, fMonomialCurrPtr, &stateNormal); // just for state normal really
	
	momentumPhase = 1.0f;    
	if ( SpinMomentum != -1 ) 
	  {
	    for ( int i = 1; i < NumberSites ; i++ )
	      {
		momentumPhase += PhaseArrayR[(i * (CurrentMomentum + SpinMomentum)) % NumberSites];
	      }
	  }
	norm = sqrtf(momentumPhase) * sqrtf(1.0f/(float)powf((float)NumberSites,(float)Filling));
			    
	int globalThreadIdx = blockIdx.x * blockDim.x + threadIdx.x;
	
	for ( int i = 0 ; i < Filling; i++ )
	  {	    
	    fMonomialCurrPtr[i] = Permutations[threadIdx.x*Filling + i];
	    if ( threadIdx.x < (blockDim.x -1) )
	      {	    
		fMonomialNextPtr[i] = Permutations[(threadIdx.x +1) * Filling + i];
	      }
	  }
		    
	int permuteFlg = 1;
	while ( permuteFlg && ((threadIdx.x == (blockDim.x -1 )) || (compare(fMonomialCurrPtr, fMonomialNextPtr, Filling) == 0)) )    
	  {
	    int oidx = blockIdx.x;    
	    while ( oidx < VectorLength ) 
	      {  
		get_monomial_fermions_cuda(BasisElements[oidx], NumberSites, Filling, oMonomialPtr);	    	    
		phaseIdx = 0;
		for ( int i = 0 ; i < Filling; i++ )
		  {    
		    phaseIdx += oMonomialPtr[i] * fMonomialCurrPtr[i];	            
		  }
		shared_tmp_vals[threadIdx.x] += PhaseArrayR[phaseIdx % NumberSites] * VectorElements[oidx] * configuration_normal_cuda(BasisElements[oidx], NumberSites, mask);
		oidx += gridDim.x;
	      }
	    permuteFlg = permute_cuda(fMonomialCurrPtr, Filling);
	  }		      
		    
	__syncthreads();
	int i = blockDim.x/2;
	while ( i != 0 ) 
	  {
	    if ( threadIdx.x < i )
	      {
		shared_tmp_vals[threadIdx.x] += shared_tmp_vals[threadIdx.x + i];
	      }
	    __syncthreads();  
	    i /= 2;
	  }
	__syncthreads();  
	if ( threadIdx.x == 0 )  
	  {
	    BlockSums[blockIdx.x] = shared_tmp_vals[0] * norm * stateNormal;
	  }
      }
  }
  
long get_rank_chunk_start(long N, int Rank, int Size){
	long  max_chunk, min_chunk, cut_off; 
	
	min_chunk = (long)floor((double)N / (double)Size); 
	max_chunk = (long)ceil((double)N / (double)Size);
	cut_off = N - (Size * min_chunk); 

	if ( Rank < cut_off ) {
		return Rank * max_chunk;
	} else {
		return (cut_off * max_chunk) + ((Rank - cut_off) * min_chunk);
	}
}  
  
void find_boundary_permutations(int *InitialSet, int Filling, int NumberSites, int NumberPartitions, int *Permutations)
  {
    //want to find how many permutations there can be. 
    int monomial[Filling];
    for ( int i = 0 ; i < Filling ; i++ )
      {
	monomial[i] = InitialSet[i];
      }
    long numberPermutations = 0; 
    
    while ( permute(monomial, Filling) ) 
      {
	numberPermutations++;
      }        
    
    //sort(monomial, monomial + Filling);    
    for ( int i = 0 ; i < Filling ; i++ )
      {
	monomial[i] = InitialSet[i];
      }
    long permutationCount = 0;
    int partitionIdx = 0;
    int lastIdx = 0;
    
    do 
      {	
	if ( permutationCount == get_rank_chunk_start(numberPermutations, partitionIdx, NumberPartitions) )
	  {
	      for ( int i = 0 ; i < Filling ; i++ ) 
		{
		  Permutations[partitionIdx * Filling + i] = monomial[i];		  
		}
	      lastIdx = partitionIdx;
	      partitionIdx++;
	  }
	permutationCount++;
      } 
    while ( permute(monomial, Filling ) && ( partitionIdx < NumberPartitions )  );    
    while ( partitionIdx < NumberPartitions ) 
      {
	for ( int i = 0 ; i < Filling ; i++ ) 
	  {
	    Permutations[partitionIdx * Filling + i] = Permutations[lastIdx * Filling + i];
	  }	
	partitionIdx++;
      }
  }


void fourier_transform_vector_cuda2(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, int CurrentMomentum, int SpinMomentum)
  {        
    //get device code
    cudaDeviceProp prop;
    int dev;
    cudaGetDevice( &dev );    
    memset( &prop, 0, sizeof( cudaDeviceProp ) );
    prop.multiProcessorCount = 16;    
    cudaChooseDevice( &dev, &prop );
    cout << "ID of CUDA device with closest to 32 cores. " << dev << endl;
    cudaSetDevice( dev );
    
    //prepare phase arrays
    float *phaseArrayR = new float[NumberSites+1];
    float *phaseArrayI = new float[NumberSites+1];
    phaseArrayR[0] = 1.0f;
    phaseArrayI[0] = 0.0f;
    for( int i = 1  ; i <= NumberSites ; i++ ) 
      {	    	    				
        float angle = ( PI * 2.0f * (float)i) / ((float)NumberSites);
	phaseArrayR[i] = (float)cos(angle);
	phaseArrayI[i] = (float)sin(angle);	
      }                           
        
    float *d_blockSums, *blockSums, *d_phaseArrayR, *d_VectorElements;
    unsigned long *d_BasisElements;
    int *d_permutations, *permutations;
    cudaError_t cudaerr;    
    permutations = new int[CUDA_THREADS_PER_BLOCK * Filling];
    
    blockSums = new float[CUDA_BLOCKS];
    cudaerr = cudaMalloc( (void**) &d_blockSums, CUDA_BLOCKS * sizeof(float)); cout << "Allocating d_blockSums: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMalloc( (void**) &d_permutations, CUDA_THREADS_PER_BLOCK * Filling * sizeof(int)); cout << "Allocating d_blockSums: " << cudaGetErrorString(cudaerr) << endl;   
    cudaerr = cudaMalloc( (void**) &d_phaseArrayR, (NumberSites + 1) * sizeof(float)); cout << "Allocating d_phaseArrayR: " << cudaGetErrorString(cudaerr) << endl;       
    //cudaerr = cudaMalloc( (void**) &d_number_permutations, CUDA_BLOCKS * CUDA_THREADS_PER_BLOCK * sizeof(int)); cout << "Allocating d_number_permutations: " << cudaGetErrorString(cudaerr) << endl;   
    cudaerr = cudaMalloc( (void**) &d_BasisElements, VectorLength * sizeof(unsigned long) );cout << "Allocating d_BasisElements: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMalloc( (void**) &d_VectorElements, VectorLength * sizeof(float) );cout << "Allocating d_VectorElements: " << cudaGetErrorString(cudaerr) << endl;
       
    cudaerr = cudaMemcpy( d_phaseArrayR, phaseArrayR, (NumberSites + 1) * sizeof(float), cudaMemcpyHostToDevice ); cout << "Copying phaseArrayR: " << cudaGetErrorString(cudaerr) << endl;    	    
    cudaerr = cudaMemcpy( d_BasisElements, BasisElements, VectorLength * sizeof(unsigned long), cudaMemcpyHostToDevice ); cout << "Copying BasisELements: " << cudaGetErrorString(cudaerr) << endl;
    cudaerr = cudaMemcpy( d_VectorElements, VectorElements, VectorLength * sizeof(float), cudaMemcpyHostToDevice ); cout << "Copying VectorElements: " << cudaGetErrorString(cudaerr) << endl;
    int monomial[Filling];    
    
    double start, finish, max, start2, finish2;
    
    max = 0.0;
    
    for ( int fidx = 0 ; fidx < FTBasisDim ; fidx++ )
    //for ( int fidx = 6 ; fidx < 7 ; fidx++ ) 
      {
	start2 = MPI_Wtime();
	get_monomial_bosons(FTBasisElements[fidx], NumberSites, Filling, monomial);
	// int count = 0;
// 	do {
// 	  cout << count++ << ": " << convert_from_monomial_bosons(NumberSites, Filling, monomial) << ": ";
// 	  for ( int i = 0 ; i < Filling ; i++ )
// 	    {
// 	      cout << monomial[i] << ", ";
// 	    }
// 	  cout << endl;
// 	} while ( permute(monomial, Filling) );
    
	FTVectorElements[fidx] = 0.0f;		
	find_boundary_permutations(monomial, Filling, NumberSites, CUDA_THREADS_PER_BLOCK, permutations); 	
	//int *number_permutations = new int[CUDA_BLOCKS * CUDA_THREADS_PER_BLOCK];	
	// for ( int i = 0 ; i < CUDA_BLOCKS * CUDA_THREADS_PER_BLOCK ; i++ ) 
// 	  {
// 	    cout << i << ": " ;
// 	    for ( int j = 0 ; j < Filling ; j++ )
// 	      {
// 		cout << permutations[i*Filling + j]  << " ";
// 	      }
// 	      cout << convert_from_monomial_bosons(NumberSites, Filling, &permutations[i*Filling]);
// 	      cout << endl;
// 	  }
	cudaerr = cudaMemcpy( d_permutations, permutations, CUDA_THREADS_PER_BLOCK * Filling * sizeof(int), cudaMemcpyHostToDevice); //cout << "Copying permutations: " << cudaGetErrorString(cudaerr) << endl;	
	
	for ( int i = 0 ; i < CUDA_BLOCKS; i++ ) 
	  {
	    blockSums[i] = 0.0f;
	  }
	cudaerr = cudaMemcpy( d_blockSums, blockSums, CUDA_BLOCKS * sizeof(float), cudaMemcpyHostToDevice); //cout << "Copying blockSums: " << cudaGetErrorString(cudaerr) << endl;  
	start = MPI_Wtime();
	calculate_FT_element<<<CUDA_BLOCKS, CUDA_THREADS_PER_BLOCK>>>(d_VectorElements, d_BasisElements, VectorLength, FTBasisElements[fidx], NumberSites, Filling, CurrentMomentum, SpinMomentum, d_phaseArrayR, d_blockSums, d_permutations);
	finish = MPI_Wtime();
	if ( (finish - start) > max) max = finish - start; 
		
	cudaerr = cudaGetLastError(); //cout << "Error after kernel execution: " << cudaGetErrorString(cudaerr) << endl;   	    
	//cudaerr = cudaMemcpy( number_permutations, d_number_permutations, CUDA_BLOCKS * CUDA_THREADS_PER_BLOCK  * sizeof(int), cudaMemcpyDeviceToHost ); //cout << "Copying back number_permutations: " << cudaGetErrorString(cudaerr) << endl;
	cudaerr = cudaMemcpy( blockSums, d_blockSums, CUDA_BLOCKS * sizeof(float), cudaMemcpyDeviceToHost );// cout << "Copying back blockSums: " << cudaGetErrorString(cudaerr) << endl;	    
	int total_permutations = 0;
	// for ( int idx = 0 ; idx < CUDA_BLOCKS ; idx++ )
// 	  {		
// 	    total_permutations = 0;
// 	    for ( int j = 0 ; j < CUDA_THREADS_PER_BLOCK ; j++ ) 
// 	      {
// 		  total_permutations += number_permutations[idx * CUDA_THREADS_PER_BLOCK + j]; 		  
// 	      }
// 	    cout << "Block: " << idx << ": " << total_permutations << endl;
// 	  }	    
// 	cout << "Fidx: " << fidx << ": " << FTVectorElements[fidx] << ", permutations: " << endl; 
	
	for ( int idx = 0 ; idx < CUDA_BLOCKS ; idx++ )
	  {
	    FTVectorElements[fidx] += blockSums[idx];
	  }
	finish2 = MPI_Wtime();
	cout << "Fidx: " << fidx << ": " << FTVectorElements[fidx] << ", time: " << finish2 - start2 << endl;  
      }        
    cout << "Max time taken for element was: " << max << endl;
        
    cudaFree(d_phaseArrayR);
    cudaFree(d_blockSums);
    cudaFree(d_permutations);
    cudaFree(d_VectorElements);    
    cudaFree(d_BasisElements);
          
    delete [] phaseArrayR;
    delete [] phaseArrayI;
    delete [] blockSums;
    delete [] permutations;
  }*/


