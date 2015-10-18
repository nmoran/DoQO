CFLAGS="-I. -Wall"

include ${SLEPC_DIR}/conf/slepc_common

all : entanglement_entropy p_entanglement_entropy p_entanglement_entropy_fourier_cut p_entanglement_entropy_batch p_entanglement_entropy_density

entanglement_entropy: entanglement_entropy.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o chkopts 
	-${CLINKER}  -o entanglement_entropy entanglement_entropy.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o ${SLEPC_LIB} 

p_entanglement_entropy: p_entanglement_entropy.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o  chkopts 
	-${CLINKER}  -o p_entanglement_entropy p_entanglement_entropy.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o  ${SLEPC_LIB} 
	
p_entanglement_entropy_fourier_cut: p_entanglement_entropy_fourier_cut.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o  chkopts 
	-${CLINKER}  -o p_entanglement_entropy_fourier_cut p_entanglement_entropy_fourier_cut.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o  ${SLEPC_LIB} 


p_entanglement_entropy_batch : p_entanglement_entropy_batch.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o  chkopts 
	-${CLINKER}  -o p_entanglement_entropy_batch p_entanglement_entropy_batch.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o  ${SLEPC_LIB} 

p_entanglement_entropy_density : p_entanglement_entropy_density.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o  chkopts 
	-${CLINKER}  -o p_entanglement_entropy_density p_entanglement_entropy_density.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o  ${SLEPC_LIB} 
