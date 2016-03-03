GFORT = gfortran
GSAC  = /Users/JFBG/src/sac/lib/sacio.a 

OBJ1  = get_spec.f
OBJ2  = pho_vars.f90 VEL_MODEL_CHECKS.f90



statsynr_gcc.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsynr_gcc.f90 $(GSAC) -o ../bin/statsynr_gcc


#
#Subroutines
#
bit_limit.x:
	$(GFORT) -c bit_limit.f90

get_spec.x:
	$(GFORT) -c get_spec.f

psplot.x:
	$(GFORT) -c psplot.f90

trace_subs.x:
	$(GFORT) -c trace_subs.f
	
transfer.x:
	$(GFORT) -c transfer.x

