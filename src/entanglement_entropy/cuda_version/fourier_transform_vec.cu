#include "host_fourier_transform.h"
#include "device_fourier_transform.h"


/*
 This function calculates the dimension of the boson basis set with given total momentum.
 
 sites : number of sites.
 filling : number of bosons.
 M : total momentum.
*/
int boson_basis_dimension(int sites, int filling, int M)
{  
  if ( sites > 0 && sites == M && filling == 1 )
    {
      return 1;
    }  
    
  if ( M < filling || M > sites*filling || sites == 0 || filling == 0 ) 
    {
      return 0; 
    } 
  else 
    {
    
      int dim = 0;  
      //use recursive procedure where can add boson and stay in position, add boson and move on or not add boson and move on.
      if ( sites <= M ) 
	{
	  dim += boson_basis_dimension(sites, filling-1, M-sites);	  
	}
      dim += boson_basis_dimension(sites-1, filling, M);
      return dim; 
    }
}


/*
 This function calculates the dimension of the boson basis set with given total momentum.
 
 sites : number of sites.
 filling : number of bosons.
 M : total momentum.
 states: array of states.
 maxm : the maximum m for a given state.
*/
int boson_basis_generate(int sites, int filling, int M, int pos, unsigned long *states, unsigned long  current_state)
{  
  if ( sites > 0 && sites == M && filling == 1 )
    {      
      current_state += 1ul << (sites + filling - 2);      
      states[pos++] = current_state;      
      return pos;
    }  
    
  if ( M < filling || M > sites*filling || sites == 0 || filling == 0 ) 
    {
      return pos; 
    } 
  else 
    {          
      //use recursive procedure where can add boson and stay in position, add boson and move on or not add boson and move on.
      pos = boson_basis_generate(sites-1, filling, M, pos, states, current_state);
      if ( sites <= M ) 
	{
	  current_state += 1ul << (sites + filling - 2);	  
	  pos = boson_basis_generate(sites, filling-1, M-sites, pos, states, current_state);
	}
      
      return pos; 
    }
}


int main ( int argc, char **argv) 
{
    int c;
    string inputFilename, basisFilename, outputPrefix; 
    int numberSites = 0, filling = 0, momentum = -1, spinMomentum = -1;
    bool truncate = false, useCuda = false, single = false, benchmark = false, benchmark2 = false;
    int rank, size;
    int cuda_device = 0;
    int benchmarkSize = 0;
    int benchmark2Size = 0;
    int offset = 0;
    
    MPI_Init (&argc, &argv);	/* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */    
     
    //Read command line arguments.
    while (1)
      {
	static struct option long_options[] =
	  {	    
	    {"vec",  required_argument,  0, 'v'},
	    {"basis",  required_argument,  0, 'b'},
	    {"num_sites",  required_argument,  0, 'n'},
	    {"filling",  required_argument,  0, 'f'},
	    {"truncate",  no_argument,  0, 't'},
	    {"momentum",  required_argument,  0, 'm'},
	    {"spin_momentum",  required_argument,  0, 'k'},
	    {"device", required_argument, 0, 'd'},
	    {"cuda",  no_argument,  0, 'c'},
	    {"single", no_argument, 0, 's'},
	    {"benchmark", required_argument, 0, 'x'},
	    {"benchmark2", required_argument, 0, 'y'},
	    {"offset", required_argument, 0, 'o'},
	    {0, 0, 0, 0}
	  };
	/* getopt_long stores the option index here. */
	int option_index = 0;
  
	c = getopt_long (argc, argv, "v:b:n:f:tm:k:csd:x:y:o:",
			long_options, &option_index);
  
	/* Detect the end of the options. */
	if (c == -1)
	  break;
  
	switch (c)
	  {	  
	    case 'v':
	      inputFilename.assign(optarg);
	      size_t found;
	      found = inputFilename.find_last_of(".");
	      outputPrefix.assign(inputFilename.substr(0,found));
	      break;  	  
	    case 'b':
	      basisFilename.assign(optarg);	      
	      break;
	    case 'n':
	      numberSites = atoi(optarg);
	      break;
	    case 'f':
	      filling = atoi(optarg);
	      break;
	    case 't':
	      truncate = true;
	      break;
	    case 'm':
	      momentum = atoi(optarg);
	      break;
	    case 'k':
	      spinMomentum = atoi(optarg);
	      break;
	    case 'c':
	      useCuda = true;
	      break;
	    case 's':
	      single = true;
	      break;
	    case 'd':
	      cuda_device = atoi(optarg);
	      break;
	    case 'x':
	      benchmark = true;
	      benchmarkSize = atoi(optarg);
	      break;  
	    case 'y':
	      benchmark2 = true;
	      benchmark2Size = atoi(optarg);
	      break;  
	    case 'o':
	      offset = atoi(optarg);
	      break;  
  
	    default:
	      MPI_Finalize();
	      abort ();
	  }
      }
      
    //Check required arguments were supplied.
    if ( inputFilename == "" ) 
      {
	cout << "No input filename supplied." << endl;
	return 1;
      }
    if ( basisFilename == "" ) 
      {
	cout << "No basis filename supplied." << endl;
	return 1;
      }
    if ( numberSites  <= 0 ) 
      {
	cout << "Invalid number of sites." << endl;
	return 1;
      }
    if ( filling  <= 0 ) 
      {
	cout << "Invalid filling." << endl;
	return 1;
      }
      
    //Print out values of arguements read.
    cout << "Input filename: " << inputFilename << endl;
    cout << "Basis filename: " << basisFilename << endl;
    cout << "Number sites: " << numberSites << endl;
    cout << "Filling: " << filling << endl;
    cout << "Truncation: " << truncate << endl;
    cout << "CUDA enabled: " << useCuda << endl;
    
    //Now read in the vector and basis. 
    //Will read in vector as just real floats for now and basis will be long integers.
    int vectorLength = 0;
    float *vectorElements;
    ifstream vectorFile(inputFilename.c_str());
    if ( vectorFile.is_open() ) 
      {
	  string line;
	  getline(vectorFile,line);
	  while ( line != ""  ) 
	    {
	      vectorLength++;
	      getline(vectorFile,line);
	    }
	  vectorFile.close();
	  vectorElements = new float[vectorLength];
	  vectorLength = 0;
	  vectorFile.open(inputFilename.c_str());
	  if ( vectorFile.is_open() ) 
	    {
		string line;
		getline(vectorFile,line);
		while ( line != ""  ) 
		  {		    		    
		    vectorElements[vectorLength++] = atof(line.c_str());
		    getline(vectorFile,line);
		  }
		vectorFile.close();
	    }
	}
    else
      {
	cout << "Problem opening " << inputFilename << endl;
	return 1;
      }
      
    int basisLength = 0;
    unsigned long *basisElements;
    ifstream basisFile(basisFilename.c_str());
    if ( basisFile.is_open() ) 
      {
	  string line;
	  getline(basisFile,line); //skip a line for the size part
	  getline(basisFile,line);
	  while ( line != "" ) 
	    {
	      basisLength++;
	      getline(basisFile,line);
	    }
	  basisFile.close();
	  basisElements = new unsigned long[basisLength];
	  basisLength = 0;
	  basisFile.open(basisFilename.c_str());
	  if ( basisFile.is_open() ) 
	    {
		string line;
		getline(basisFile,line); //skip a line for the size part
		getline(basisFile,line);
		while ( line != ""  ) 
		  {		    		    
		    basisElements[basisLength++] = atol(line.c_str());
		    getline(basisFile,line);
		  }
		basisFile.close();
	    }
	}
    else
      {
	cout << "Problem opening " << inputFilename << endl;
	return 1;
      }
      
    if ( vectorLength != basisLength ) 
      {
	cout << "Vector and basis lengths do not match." << endl;
	return 1;
      }
    
    cout << "Input vector size: " << vectorLength << endl;
      
    // cout << "Vector and basis details." << endl;
//     for ( int i = 0 ; i < vectorLength ; i++ ) 
//       {
// 	cout << basisElements[i] << "\t: " << vectorElements[i] << endl;
//       }
      
    int startMomentum = 0; 
    int endMomentum = numberSites * filling;
    if ( momentum >= 0 && momentum <= (numberSites * filling) ) 
      {
	startMomentum = momentum;
	endMomentum = momentum;
      }
      
    //if we are using translational invariance then apply normal now.  
    if ( spinMomentum != -1 ) 
      {
	unsigned long mask = (1ul << numberSites) - 1ul;
	for ( int i = 0 ; i < vectorLength ; i++ ) 
	  {
	    vectorElements[i] *= configuration_normal(basisElements[i], numberSites, mask);
	  }  
      }
    
    for ( int currentMomentum = startMomentum ; currentMomentum <= endMomentum ; currentMomentum++ )
      {
	//Now need to calculate the basis that the fourier transformed vector will use.
	int FTBasisDim = boson_basis_dimension(numberSites, filling, currentMomentum);
	unsigned long * FTBasisElements = new unsigned long [FTBasisDim];
	float* FTVectorElements = new float[FTBasisDim];
	for ( int i = 0 ; i < FTBasisDim ; i++ ) 
	  {
	    FTVectorElements[i] = 0.0;
	  }
	boson_basis_generate(numberSites, filling, currentMomentum, 0, FTBasisElements, 0ul);
	cout << "Sector with total momentum " << currentMomentum << " has " << FTBasisDim << " elements." << endl;
		
	if ( benchmark ) 
	  {
	    if ( benchmarkSize > 0 ) 
	      {
		unsigned long *FTBasisBenchmarkElements = new unsigned long [benchmarkSize];
		int bidx = 0, fidx = 0;
		float normal;
		int monomial[filling];
		while ( (bidx < benchmarkSize) && (fidx < FTBasisDim) )
		  {
		    get_monomial_bosons(FTBasisElements[fidx], numberSites, filling, monomial, &normal);
		    if ( normal == 1.0 ) 
		      {
			FTBasisBenchmarkElements[bidx++] = FTBasisElements[fidx];
		      }
		    fidx++;
		  }	
		if ( bidx < benchmarkSize ) 
		  {
		    cout << "Only using " << bidx << " elements for benchmark out of requested " << benchmarkSize << endl;
		    benchmarkSize = bidx; 
		  }
		if ( useCuda ) 
		  {
		    double start, finish;	   
		    start = MPI_Wtime();
		    if ( benchmark2 && benchmark2Size <= vectorLength ) 
		      {
			fourier_transform_vector_cuda(benchmark2Size, vectorElements, basisElements, benchmarkSize, FTVectorElements, FTBasisBenchmarkElements, numberSites, filling, currentMomentum, spinMomentum, single, cuda_device, outputPrefix, offset);
		      }
		    else
		      {
			fourier_transform_vector_cuda(vectorLength, vectorElements, basisElements, benchmarkSize, FTVectorElements, FTBasisBenchmarkElements, numberSites, filling, currentMomentum, spinMomentum, single, cuda_device, outputPrefix, offset);	
		      }
		    finish = MPI_Wtime();  		
		    cout << "Fourier transform of " << benchmarkSize << " benchmark elements using CUDA for momentum " << currentMomentum << " took " << finish - start << " seconds." << endl;
		  } 
		else 
		  {
		    double start, finish;	   
		    start = MPI_Wtime();
		    if ( benchmark2 && benchmark2Size <= vectorLength ) 
		      {
			fourier_transform_vector(benchmark2Size, vectorElements, basisElements, benchmarkSize, FTVectorElements, FTBasisBenchmarkElements, numberSites, filling, currentMomentum, spinMomentum);
		      }
		    else
		      {
			fourier_transform_vector(vectorLength, vectorElements, basisElements, benchmarkSize, FTVectorElements, FTBasisBenchmarkElements, numberSites, filling, currentMomentum, spinMomentum);	
		      }
		    finish = MPI_Wtime();  		
		    cout << "Fourier transform of " << benchmarkSize << " benchmark elements for momentum " << currentMomentum << " took " << finish - start << " seconds." << endl;
		  }
		delete []FTBasisBenchmarkElements;
	      }	    
	  }
	else 
	  {
	    if ( useCuda ) 
	      {
		double start, finish;	   
		start = MPI_Wtime();
		fourier_transform_vector_cuda(vectorLength, vectorElements, basisElements, FTBasisDim, FTVectorElements, FTBasisElements, numberSites, filling, currentMomentum, spinMomentum, single, cuda_device, outputPrefix, offset);	
		finish = MPI_Wtime();
		cout << "Fourier transform using CUDA for momentum " << currentMomentum << " took " << finish - start << " seconds." << endl;
	      }
	    else
	      {				
		double start, finish;
		start = MPI_Wtime();
		fourier_transform_vector(vectorLength, vectorElements, basisElements, FTBasisDim, FTVectorElements, FTBasisElements, numberSites, filling, currentMomentum, spinMomentum);	
		finish = MPI_Wtime();
		cout << "Fourier transform for momentum " << currentMomentum << " took " << finish - start << " seconds." << endl;
		
	      }
	  }	      	      	
	
	float mag = 0;
	for( int i = 0 ; i < FTBasisDim ; i++ )
	  {
	    mag += FTVectorElements[i] * FTVectorElements[i];
	  }
	cout.precision(14);
	cout << scientific << "FT vector magnitude: " << sqrt(mag) << endl;
// 	for( int i = 0 ; i < FTBasisDim ; i++ )
// 	  {
// 	    cout.precision(14);
// 	    cout << scientific << "Element " << i << ": " << FTVectorElements[i] << endl;
// 	  }  	

	//output basis and vector elements. 
	string outputFilename;
	stringstream ss("");
	ss << outputPrefix << "_FT_M_" << currentMomentum << ".vec" ; 
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
	ss.str("");
	ss << outputPrefix << "_FT_M_" << currentMomentum << ".basis" ; 
	out.open(ss.str().c_str());
	if ( out.is_open() )
	  {	    
	    for ( int i = FTBasisDim - 1 ; i >= 0 ; i-- ) 
	      {
		out << FTBasisElements[i] << endl;
	      }
	    out.close();
	  }	  	
	delete [] FTBasisElements;
	delete [] FTVectorElements;
      }  
      
    delete [] basisElements;
    delete [] vectorElements;
    MPI_Finalize();
  
    return 0;
}