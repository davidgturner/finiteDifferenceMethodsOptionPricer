all:
	g++ FDM_main.cpp ImplicitFDM.cpp ExplicitFDM.cpp CN_FDM.cpp \
	ImplicitSORFDM.cpp CN_SORFDM.cpp FDM_utils.cpp \
	BlackScholesFormula.cpp -o FDM
