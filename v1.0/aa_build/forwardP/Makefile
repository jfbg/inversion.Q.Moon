GFORT = gfortran
GSAC  = /usr/local/sac/lib/sacio.a 

OBJ1  = get_spec.f
OBJ2  = pho_vars.f90 VEL_MODEL_CHECKS.f90
OBJ3  = pho_vars.f90 VEL_MODEL_CHECKS_crfl.f90


statsynr_gcc.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsynr_gcc.f90 $(GSAC) -o ../bin/statsynr_gcc

statsynr_gcc_FirstS.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsynr_gcc_FirstS.f90 $(GSAC) -o ../bin/statsynr_gcc_FirstS

statsynr_gcc_BM_Earth.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ3) statsynr_gcc_BM_EARTH.f90 $(GSAC) -o ../bin/statsynr_gcc_BM_Earth

statsynr_gcc_BM_Moon.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ3) statsynr_gcc_BM_MOON.f90 $(GSAC) -o ../bin/statsynr_gcc_BM_Moon

statsynr_BASIN_gcc.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsynr_BASIN_gcc.f90 $(GSAC) -o ../bin/statsynr_BASIN_gcc

statsyn_gcc_TRACKPARAM.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_gcc_TRACKPARAM.f90 $(GSAC) -o ../bin/statsyn_gcc_TRACKPARAM

statsyn_gcc_COUNTHITS.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_gcc_COUNTHITS.f90 $(GSAC) -o ../bin/statsyn_gcc_COUNTHITS

statsyn_gcc_OUTPUTANGLE.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_gcc_OUTPUTANGLE.f90 $(GSAC) -o ../bin/statsyn_gcc_OUTPUTANGLE

generatespikes.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) generatespikes.f90 $(GSAC) -o ../bin/generatespikes

statsyn_gcc_OnlySArrival.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_gcc_OnlySArrival.f90 $(GSAC) -o ../bin/statsyn_gcc_OnlySArrival

statsyn_gcc_2dscat.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_gcc_2dscat.f90 $(GSAC) -o ../bin/statsyn_gcc_2dscat

statsyn_gcc_track.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_gcc_track.f90 $(GSAC) -o ../bin/statsyn_gcc_track

statsyn_global.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat.f90 $(GSAC) -o ../bin/statsyn_globalscat

statsyn_track.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_track_gcc.f90 $(GSAC) -o ../bin/statsyn_track

statsyn_global_CRFL_noATT.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_CRFL_noATT.f90 $(GSAC) -o ../bin/statsyn_globalscat_CRFL_noATT

statsyn_global_CRFL.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_CRFL.f90 $(GSAC) -o ../bin/statsyn_globalscat_CRFL

statsyn_global_BM_AMP.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_BM_AMP.f90 $(GSAC) -o ../bin/statsyn_globalscat_BM_AMP

statsyn_global_BM.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_BM.f90 $(GSAC) -o ../bin/statsyn_globalscat_BM

trackglobal.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat.f90 $(GSAC) -o ../bin/statsyn_global
    
trackglobal_Energy.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_Energy.f90 $(GSAC) -o ../bin/statsyn_global
    
trackglobal_FlatRho.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_FlatRho.f90 $(GSAC) -o ../bin/statsyn_global

statsyn_VPREM.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_VPREM.f90 $(GSAC) -o ../bin/statsyn_VPREM
    
statsyn_VPREM_BEN.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_VPREM_BEN.f90 $(GSAC) -o ../bin/statsyn_VPREM_BEN

statsyn_VPREM_BEN_vs.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_VPREM_BEN_vs.f90 $(GSAC) -o ../bin/statsyn_VPREM_BEN_vs

statsyn_VPREM_NORM.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_VPREM_NORM.f90 $(GSAC) -o ../bin/statsyn_VPREM_NORM

trackglobal_FlatRho_isoE.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_FlatRho_isoE.f90 $(GSAC) -o ../bin/statsyn_global_isoE
	
trackglobal_FlatRho_isoE_POWER.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_FlatRho_isoE_Power.f90 $(GSAC) -o ../bin/statsyn_global_isoE_POWER
    
trackglobal_OnlyP.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_OnlyP.f90 $(GSAC) -o ../bin/statsyn_global

trackglobalsingle.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_globalscat_single.f90 $(GSAC) -o ../bin/statsyn_global_single

trackiso.x:
	$(GFORT) -m64 $(OBJ1) $(OBJ2) statsyn_TRACK_iso.f90 $(GSAC) -o ../bin/statsyn_iso

tracktime.x:	
	$(GFORT) -m64 get_spec.f statsyn_1D_TIMETRACK.f90 $(GSAC) -o ../bin/statsyn_1D_track_time

track.x:	
	$(GFORT) -m64 get_spec.f statsyn_TRACK.f90 $(GSAC) -o ../bin/statsyn_TRACK

statsyn_1D.x:
	$(GFORT) -m64 get_spec.f statsyn_1D.f90 $(GSAC) -o ../bin/statsyn_1D

statsyn_1D_att.x:
	$(GFORT) get_spec.f statsyn_1D_att.f90 $(GSAC) -o ../bin/statsyn_1D_att

trace.nw.o:
	$(FORT) trace.nw.f -o ../bin/trace.nw.1

get_spec.o:
	$(FORT) -c get_spec.f

pstack2D.x:
	$(GFORT) pstack2D.f90  psplot.f90 -o ../bin/pstack2D

interp_mod.o:
	$(GFORT) interp_mod.f -o ../bin/interp_mod

decon.o:
	$(GFORT) decon.f transfer.f get_spec.f -o ../bin/decon

add_noise.x:
	$(GFORT) add_noise.f90 transfer.f get_spec.f bit_limit.o -o add_noise

convolve.o:
	$(GFORT) convolve.f transfer.f get_spec.f bit_limit.o -o ../bin/convolve

trace.luna.o:
	fort77 trace.luna.f trace_subs.f psplot.f $(XLIB) -o ../bin/trace.luna

xyz2sac.x:
	$(GFORT) xyz2sac.f90 $(GSAC) -o ../bin/xyz2sac

sac2xyz.x:
	$(GFORT) sac2xyz.f90 $(GSAC) -o ../bin/sac2xyz

convert_prem_moon.x:
	$(GFORT) convert_prem_moon.f90 -o ../bin/convert_prem_moon


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

