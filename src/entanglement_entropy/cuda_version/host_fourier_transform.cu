#include "host_fourier_transform.h"

#define PI 3.141592653589793f

// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
void get_monomial_bosons(unsigned long Index, int Sites, int Bosons, int *Monomial) 
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
void get_monomial_bosons(unsigned long Index, int Sites, int Bosons, int *Monomial, float *Normal) 
  {
    int mom = 1;
    int idx = 0;
    int pos = 0;
    *Normal = 1.0;
    double num = 0;
    while ( pos < (Sites + Bosons) ) {
      if ( (Index & (1ul << pos) ) > 0ul ) {
	  Monomial[idx++] = mom;
	  num += 1.0;
	  *Normal *= num;
      } else {
	  mom += 1;
	  //if ( num > 1.0 ) *normal *= num; 
	  num = 0.0; 
      }	  
      pos++;
    }
    //if ( num > 1.0 ) *normal *= num; 
    *Normal = sqrt(*Normal);
  }


// calculated the momentum that each boson is at and populates the monomials array with this. Does so in ascending order.
void get_monomial_fermions(unsigned long Index, int Sites, int Fermions, int *Monomial) 
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

int permute(int *Monomial, int Len)
  {
    int k,l,tmp;
    
    k = Len - 2;
    while (k > 0 && Monomial[k] >= Monomial[k+1] ) k--;
    
    if ( k == 0 && (Monomial[0] >= Monomial[1]) ) return 0;
    
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

float configuration_normal(unsigned long Config, int NumberSites, unsigned long Mask)
  {
    float normal = 1.0;    
    
    for ( int i = 1 ; i < NumberSites ; i++ ) 
    {
      if ( Config == ( ((Config << i) | (Config >> (NumberSites - i))) & Mask ) ) 
	{
	  normal += 1.0;
	}
    }
    return 1.0/sqrt(normal);             
  }

void fourier_transform_vector(int VectorLength, float *VectorElements, unsigned long *BasisElements, int FTBasisDim, float *FTVectorElements, unsigned long *FTBasisElements, int NumberSites, int Filling, int CurrentMomentum, int SpinMomentum)
  {
    int *fmonomial, *fmonomialTmp, *omonomial;
    fmonomial = new int[Filling];
    fmonomialTmp = new int[Filling];
    omonomial = new int[Filling];
    
    float *phaseArrayR = new float[NumberSites+1];
    float *phaseArrayI = new float[NumberSites+1];
    phaseArrayR[0] = 1.0;
    phaseArrayI[0] = 0.0;
    for( int i = 1  ; i <= NumberSites ; i++ ) 
      {	    	    				
        double angle = ( PI * 2.0 * (double)i) / ((double)NumberSites);
	phaseArrayR[i] = cos(angle);
	phaseArrayI[i] = sin(angle);
	//cout << "Phase array[" << i << "] = " << phaseArrayR[i] << endl;
      }
    
    float norm = sqrt(1.0/pow((double)NumberSites,(double)Filling));
    
    float momentumPhase = 1.0;
    if ( SpinMomentum != -1 ) 
      {
	for ( int i = 1; i < NumberSites ; i++ )
	  {
	    momentumPhase += phaseArrayR[(i * (CurrentMomentum + SpinMomentum)) % NumberSites];
	  }
      }
    norm *= sqrt(momentumPhase);
    
    int phaseIdx;
    float tmp_val, tmp_val2;
    float stateNormal = 1.0;
    
    if ( SpinMomentum == -1 ) 
      {
	for ( int fidx = 0 ; fidx < FTBasisDim ; fidx++ )
	  {
	    tmp_val = 0.0;
	    get_monomial_bosons(FTBasisElements[fidx], NumberSites, Filling, fmonomial, &stateNormal);
	    get_monomial_bosons(FTBasisElements[fidx], NumberSites, Filling, fmonomialTmp, &stateNormal);
	    for ( int oidx = 0 ; oidx < VectorLength ; oidx++ ) 
	      {
		//sort(fmonomial, fmonomial + Filling); //reset to original configuration.
		memcpy(fmonomial, fmonomialTmp, sizeof(int)*Filling);
		tmp_val2 = 0.0;
		get_monomial_fermions(BasisElements[oidx], NumberSites, Filling, omonomial);
		do 
		  {
		    phaseIdx = 0;
		    for ( int i = 0 ; i < Filling; i++ )
		      {
			phaseIdx += omonomial[i] * fmonomial[i];	            
		      }
		    tmp_val2 += phaseArrayR[phaseIdx % NumberSites];
		  }
		while ( permute(fmonomial, Filling) );
		tmp_val += tmp_val2 * VectorElements[oidx];
	      }	
	    FTVectorElements[fidx] = tmp_val * norm * stateNormal;
	  }
      }
    else 
      {
	unsigned long mask = (1ul << NumberSites) - 1ul;
	for ( int fidx = 0 ; fidx < FTBasisDim ; fidx++ )
	  {
	    tmp_val = 0.0;
	    get_monomial_bosons(FTBasisElements[fidx], NumberSites, Filling, fmonomial, &stateNormal);
	    get_monomial_bosons(FTBasisElements[fidx], NumberSites, Filling, fmonomialTmp, &stateNormal);
	    for ( int oidx = 0 ; oidx < VectorLength ; oidx++ ) 
	      {
		//sort(fmonomial, fmonomial + Filling); //reset to original configuration.
		memcpy(fmonomial, fmonomialTmp, sizeof(int)*Filling);
		tmp_val2 = 0.0;
		get_monomial_fermions(BasisElements[oidx], NumberSites, Filling, omonomial);
		do 
		  {
		    phaseIdx = 0;
		    for ( int i = 0 ; i < Filling; i++ )
		      {
			phaseIdx += omonomial[i] * fmonomial[i];	            
		      }
		    tmp_val2 += phaseArrayR[phaseIdx % NumberSites];
		  }
		while ( permute(fmonomial, Filling) );
		tmp_val += tmp_val2 * VectorElements[oidx] * configuration_normal(BasisElements[oidx], NumberSites, mask);
	      }	
	    FTVectorElements[fidx] = tmp_val * norm * stateNormal;
	  }
      }
    delete [] fmonomial;
    delete [] fmonomialTmp;
    delete [] omonomial;
    delete [] phaseArrayR;
    delete [] phaseArrayI;
  }
