      subroutine forwardP(vp, vs, zs, nlay, nbs, distances, &
     &   T1s, T2s, h1s, h2s,  traveltimes,vpc,vpm,vsc,vsm)
!subroutine getTTsurface
! ---------------------------------------------------------------------
!  DECLARATIONS

!   nlay = nombre de couches dans le modèle
!        nbs = nombre d'observations
!   distances = distances épicentrales (nbs)

      use reloc_vars   !global variables
      
      implicit none
      
!      real       T1s(1000),T2s(1000),h1s(1000), h2s(1000)
!      real       distances(*), traveltimes(1000)
      integer  :: nbs  !number of distances to compute time for.
!      real(8)  T1s(*), T2s(*), h1s(*), h2s(*), distances(*), traveltimes(*)
      real(8)  T1s(*), T2s(*), h1s(*), h2s(*)
      real(8)    traveltimes(*),distances(*)
      real       profmin, T1, h1, T2, h2
      real       vpm, vpc, vsm, vsc
      integer    length_thin, length_thick, II, JJ, DD
      real       sourcedepth
      
      real       thin_depth(6), thin_vp(6), thin_vs(6)
      real       thick_depth(6), thick_vp(6), thick_vs(6)
      
     
      integer    numdist, pathTYPE, pathTYPE_INIT
      

	  real(8)      PAlldistB(nump0,maxB), PAllttB(nump0,maxB), PAllraypB(nump0,maxB)
      real(8)      SAlldistB(nump0,maxB), SAllttB(nump0,maxB), SAllraypB(nump0,maxB)  
      integer      PlengthB(maxB), PbitB    
      integer      SlengthB(maxB), SbitB    
      
            
      !SURFACETT
      real       maxdist,deg2km
      real   outdist(181)
      real   outttvp_thin(181), outrayp_thin(181)
      real   outttvs_thin(181), outrays_thin(181)
      real   outttvp_thick(181), outrayp_thick(181)
      real   outttvs_thick(181), outrays_thick(181)
      real   thin_cm, thick_cm !depth of crust-mantle boundary if surface is 0km
      
      real    TTsum(181) 
      real(8) TTdist,tempdist,TTdistPARTIAL
      real   epidist
      
      real   disttest
      
      integer check_shadows
  
  
	  !INTERIORTT

!       integer, parameter :: npick0= 500, nq0=1000  !also in locate_events
	  integer, parameter :: raypnum = 400
	  real    outttINT(raypnum),outraypINT(raypnum),outdistINT(raypnum)
	  real    outraypINTcorr(raypnum)

	  real(8)    INTdist(raypnum),INTtt(raypnum)
	  real(8)    H1dist(raypnum),H1tt(raypnum)
	   real(8)    H2dist(raypnum),H2tt(raypnum)
	  real(8)    PARTIALdist(raypnum),PARTIALtt(raypnum)
	  real    raypmin,raypmax,raypstep
	  integer index,raypnum2, fullinterior




	  real    ALLrayp(raypnum)
	  real(8) ALLdist(raypnum),ALLtt(raypnum)
	  integer  indexB,dirB,bitB, flagPASS
	  real(8)               AlltempB(raypnum)
	  integer  lengthB(150)
	  real(8) FINALrayp, FINALdist, FINALtt,FINALraypB(150),FINALttB(150)
	  real(8) FINALtt_SURF

			   !surfaceINT
	  real    outttH1(raypnum),outraypH1(raypnum),outdistH1(raypnum),outraypH1corr(raypnum)
	   real    outttH2(raypnum),outraypH2(raypnum),outdistH2(raypnum),outraypH2corr(raypnum)
	   real    outttPARTIAL(raypnum),outraypPARTIAL(raypnum),outdistPARTIAL(raypnum)
 
	  !TEMP
	  integer, parameter ::  nlay0=4000
	  integer             status, nlay, nbMo, i
	  !TEMP


	  real(8)        vp(*),vs(*),zs(*)          

	  real(8)                     checkNAN
	  integer                     check99   
	  integer           mindistindex, KK, j
	  real(8)           mindistval   
	  
	  real         Pdel123(3), Pdepth123(3)
	  real(8)      Pptab8(nump0), Ptttab8(nump0), Pdeltab8(nump0),Pztab8(nump0)
	  integer      Pndel,Pndepth
	  real         Sdel123(3), Sdepth123(3)
	  real(8)      Sptab8(nump0), Stttab8(nump0), Sdeltab8(nump0),Sztab8(nump0)
	  integer      Sndel,Sndepth	  
      real, allocatable, dimension(:,:)  ::  Ptt, Stt
      real, allocatable, dimension(:)    ::  deltab, deptab
!       character*250 file_picks
!       real qlat2(nq0),qlon2(nq0),qdep2(nq0),qorg2(nq0)
!       integer    ipick
      real       deltaf(npick0)
      real       resid(npick0)
!       integer*2  idph(npick0),idsta(npick0)
      real rms(nq0),rmed(nq0),resol(nq0)
!       integer   nq
!       real       rms_all
      
      real pi
      parameter (pi=3.1415927)


      ecircum=2.*pi*radius
      kmdeg=ecircum/360.
      degrad=180./pi

! ---------------------------------------------------------------------


! COMPUTE P-WAVE TRAVEL TIMES
        Pdel123(1) = 0
        Pdel123(2) = 180
        Pdel123(3) = 2
        Pndel = floor((Pdel123(2)-Pdel123(1))/Pdel123(3))+1
        Pdepth123(1) = 0
        Pdepth123(2) = 1
        Pdepth123(3) = 1
        Pndepth = floor((Pdepth123(2)-Pdepth123(1))/Pdepth123(3))+1  
        allocate( Ptt(Pndel,Pndepth) ) 
        allocate( deptab(Pndepth) ) 
        allocate( deltab(Pndel) )                         
        do KK = 1,Pndel
         deltab(KK) = (KK-1) * Pdel123(3) + Pdel123(1)
        end do
        do KK = 1,Pndepth
         deptab(KK) = (KK-1) * Pdepth123(3) + Pdepth123(1)
        end do
        
! COMPUTE S-WAVE TRAVEL TIMES        
        Sdel123(1) = 0
        Sdel123(2) = 180
        Sdel123(3) = 2
        Sndel = floor((Sdel123(2)-Sdel123(1))/Sdel123(3))+1
        Sdepth123(1) = 0
        Sdepth123(2) = 1
        Sdepth123(3) = 1
        Sndepth = floor((Sdepth123(2)-Sdepth123(1))/Sdepth123(3))+1   
        allocate( Stt(Sndel,Sndepth) ) 


		call compute_ttimes(1, vp, vs, zs, nlay, &
     &Pptab8, Ptttab8, Pdeltab8, Pztab8, Ptt, Pndel, Pdel123, Pndepth, Pdepth123)
     
     !Pdeltab8 in k, convert to deg
        Pdeltab8 = Pdeltab8 / kmdeg   
        

    
		call compute_ttimes(2, vp, vs, zs, nlay, &
     &Sptab8, Stttab8, Sdeltab8, Sztab8, Stt, Sndel, Sdel123, Sndepth, Sdepth123)     

     !Pdeltab8 in k, convert to deg
        Sdeltab8 = Sdeltab8 / kmdeg   


        do ii = 1,nump
          write(301,*)  Pptab8(ii), Ptttab8(ii), Pdeltab8(ii), Pztab8(ii)
          write(302,*)  Sptab8(ii), Stttab8(ii), Sdeltab8(ii), Sztab8(ii)        
        end do
        
        
!         do ii = 1,nump
!         print*,Pdeltab8(ii),Sdeltab8(ii),Pptab8(ii),Sptab8(ii),Ptttab8(ii),Stttab8(ii)
!         end do

      call locate_events(Ptt,Stt,Pndel,deltab,Pndepth,&
     &deptab,rms,&
     &rmed,resol,deltaf,resid)
     
!      do ii = 1,ipick
!      print*,'H',ii,deltaf(ii)
!      end do

     if (do_shadows == 1) then

		  call break_ttcurve(1,vp,Pdeltab8,Ptttab8,Pptab8,PAlldistB,&
		 &PAllttB,PAllraypB,PlengthB,PbitB,traveltimes,deltaf)
	  
		  call break_ttcurve(2,vs,Sdeltab8,Stttab8,Sptab8,SAlldistB,&
		 &SAllttB,SAllraypB,SlengthB,SbitB,traveltimes,deltaf)     

     end if

	! Compute rms of all locations     
     rms_all = 0.
     check_shadows = 0    
     do kk = 1,ipick
        rms_all = rms_all + (resid(kk))**2
        if (do_shadows == 1) then
	        distances(kk) = real(deltaf(kk),8)
    	    if (traveltimes(kk) > 1e8 ) check_shadows = 1
    	endif
     end do
     rms_all = (rms_all/real(ipick))**.5
     
     if (do_shadows == 1) then
	     if (check_shadows == 1) rms_all = 1e9
	 end if
	 
	 if (build_synth == 1) then
	     call write_synth_data(vp,vs,Pdeltab8,Ptttab8,Pptab8,Sdeltab8,Stttab8,Sptab8)
	 end if
     
   
    
       end subroutine forwardP



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE break_ttcurve(phase,vel,Alldist2,Alltt2,Allrayp2,AlldistB,AllttB,AllraypB,&
     &lengthB,bitB,traveltimes,deltaf)


       use reloc_vars
       
       implicit none

       
	   real(8)      ALLdist2(nump0),ALLtt2(nump0), ALLrayp2(nump0)
	   integer      indexB,dirB,bitB, flagPASS
	   real(8)      AlldistB(nump0,maxB), AllttB(nump0,maxB), AllraypB(nump0,maxB)
	   real(8)      AlltempB(nump0)
	   integer      lengthB(maxB)
	   real(8)      checkNAN
	   integer      check99   
	   integer      mindistindex, KK, minraypindex
	   real(8)      mindistval  
       real(8)      FINALrayp, FINALdist, FINALtt,FINALraypB(maxB),FINALttB(maxB)
       integer      II, PP, WW, QQ
       real(8)      vel(*)
       integer      phase
       real(8)      traveltimes(*)
       real         epidist,deltaf(*)

     	   

       if ((phase.ne.1).and.(phase.ne.2)) THEN
    	   print*, 'Phase must be 1 (p) or 2 (S), now is:',phase
           stop
       end if
       

	   dirB = int((ALLdist2(2)-Alldist2(1))/abs(ALLdist2(2)-Alldist2(1)))
	   bitB = 1
	   indexB = 3

	   AlldistB = -9999
	   AllttB = -9999
	   AllraypB = -9999

	   AlldistB(1,1) = Alldist2(1)
	   AllttB(1,1) = Alltt2(1)
	   AllraypB(1,1) = Allrayp2(1)
	   AlldistB(2,1) = Alldist2(2)
	   AllttB(2,1) = Alltt2(2)
	   AllraypB(2,1) = Allrayp2(2)

	   DO II = 3,nump-1
	   
! 	     print*, Alldist2(II), alltt2(II), allrayp2(II), 1/vel(1)


		 IF (int((ALLdist2(II)-Alldist2(II-1))/abs(ALLdist2(II)-Alldist2(II-1))).ne.dirB) THEN

		   IF (dirB.eq.-1) THEN  !Reverse bitB to get increasing distances

			 AlldistB(1:indexB-3,bitB) = AlldistB(indexB-2:2:-1,bitB)
			 AllttB(1:indexB-3,bitB) = AllttB(indexB-2:2:-1,bitB)
			 AllraypB(1:indexB-3,bitB) = AllraypB(indexB-2:2:-1,bitB)

		   ELSE

			 AlldistB(1:indexB-3,bitB) = AlldistB(2:indexB-2,bitB)
			 AllttB(1:indexB-3,bitB) = AllttB(2:indexB-2,bitB)
			 AllraypB(1:indexB-3,bitB) = AllraypB(2:indexB-2,bitB)

		   END IF

		   lengthB(bitB) = indexB-3

		   checkNAN = 0.
		   check99 = 0
		   DO KK = 1,lengthB(bitB)
			 checkNAN = checkNAN + AlldistB(KK,bitB)
			 if (AlldistB(KK,bitB).eq.-9999) THEN
			   check99 = 1;
			 end if   
		   END DO

		   IF ((isnan(checkNAN)).or.(check99.eq.1)) THEN
		   ELSE     
			 IF (checkNAN.ne.0) THEN
			   bitB = bitB + 1
			 END IF         
		   END IF


		   dirB = int((ALLdist2(II)-Alldist2(II-1))/abs(ALLdist2(II)-Alldist2(II-1)))
		   indexB = 3

		   AlldistB(1,bitB) = Alldist2(II-1)
		   AllttB(1,bitB) = Alltt2(II-1)
		   AllraypB(1,bitB) = Allrayp2(II-1)
		   AlldistB(2,bitB) = Alldist2(II)
		   AllttB(2,bitB) = Alltt2(II)
		   AllraypB(2,bitB) = Allrayp2(II)
				   
		 ELSE

		   AlldistB(indexB,bitB) = Alldist2(II)
		   AllttB(indexB,bitB) = Alltt2(II)
		   AllraypB(indexB,bitB) = Allrayp2(II)

		   indexB = indexB + 1


		 END IF


	   END DO


	   IF (dirB.eq.-1) THEN  !Reverse bitB to get increasing distances

		 AlldistB(1:indexB-3,bitB) = AlldistB(indexB-2:2:-1,bitB)
		 AllttB(1:indexB-3,bitB) = AllttB(indexB-2:2:-1,bitB)
		 AllraypB(1:indexB-3,bitB) = AllraypB(indexB-2:2:-1,bitB)

	   END IF

	   lengthB(bitB) = indexB-3

	   IF (lengthB(bitB).eq.0) THEN   !Check if length is 0
		 bitB = bitB -1
	   END IF

	   IF (lengthB(bitB).eq.2) THEN   !Check if length is 2 and two are the same
		 IF ((AlldistB(1,bitB) + AlldistB(2,bitB)).eq.0) THEN
		   bitB = bitB -1
		 END IF
	   END IF
	   
! 	   print*,bitB,lengthB(bitB),'HAHA'

	   mindistindex = 1
	   mindistval = AlldistB(1,1)

	   !!  FIX short distances 0deg
	   IF (vel(2).gt.vel(1)) THEN     !Only do if gradient is positif at surface
		 DO II = 2,bitB
		   IF (AlldistB(1,II).lt.AlldistB(1,mindistindex)) THEN
			 mindistindex = II
			 mindistval = AlldistB(1,II)
		   END IF
		 END DO

		 !Distance   
		 AlltempB = AlldistB(:,mindistindex)
		 AlldistB(1,mindistindex) = 0.
		 AlldistB(2:lengthB(mindistindex)+1,mindistindex) = AlltempB(1:lengthB(mindistindex))
		 !TT
		 AlltempB = AllttB(:,mindistindex)
		 AllttB(1,mindistindex) = 0.
		 AllttB(2:lengthB(mindistindex)+1,mindistindex) = AlltempB(1:lengthB(mindistindex))
		 !rayp
		 AlltempB = AllraypB(:,mindistindex)
		 AllraypB(1,mindistindex) = 1/vel(1)     !pmax
		 AllraypB(2:lengthB(mindistindex)+1,mindistindex) = AlltempB(1:lengthB(mindistindex))   

		 lengthB(mindistindex) = lengthB(mindistindex) + 1

	   END IF
	   
	   
	   
	   DO PP = 1,ipick
	   
	   FINALraypB = 1e9
	   FINALttB = 1e9	   
	   
	       IF (idph(PP).eq.phase) THEN
	       
           epidist = deltaf(PP)
	   
		   DO II = 1,bitB
		  
	!          print*,II,ALLdistB(1,II),ALLdistB(lengthB(II),II)
		  
			  IF ((epidist.ge.ALLdistB(1,II)).and.(epidist.le.ALLdistB(lengthB(II),II))) THEN
	!          print*,'JAJAJ'
				call interp_linear(1,lengthB(II),ALLdistB(1:lengthB(II),II),ALLraypB(1:lengthB(II),II), &
							  1,real(epidist,8),FINALraypB(II))
				call interp_linear(1,lengthB(II),ALLdistB(1:lengthB(II),II),ALLttB(1:lengthB(II),II),  &
							  1,real(epidist,8),FINALttB(II))
			  END IF
			
	!          print*,II,FINALraypB(II),FINALttB(II)
			END DO

!              if (PP.eq.79) then
!              	   DO QQ = 1,bitB
! 				 print*, lengthB(QQ)
! 				   do WW = 1,lengthB(QQ)
! 					write(100+PP,*) WW, ALLdistB(WW,QQ), ALLttB(WW,QQ)
! 				 end do		
! 				end do
!                 print*,'AHAHAHAHA',FINALttB(1:bitB), epidist
!              end if
              

			  IF (interpMETHOD.eq.1) THEN !minimize rayp
				  minraypINDEX =  minloc(FINALraypB(1:bitB),1)
				  FINALtt = FINALttB(minraypINDEX)
			  ELSE IF (interpMETHOD.eq.2) THEN !minimize tt
				  FINALtt = minval(FINALttB(1:bitB))
			  END IF
		  

		 
			 traveltimes(PP) = FINALtt
			 
! 			  write(*,*)   PP,deltaf(PP),idph(PP),traveltimes(PP),qindex(PP)
	    
	         END IF
       END DO
	   
	   
	   
!        DO II = 1,bitB
!          print*,'bit:',II,lengthB(II)
!        END DO	   


END SUBROUTINE break_ttcurve


SUBROUTINE break_ttcurve_synth(phase,vel,Alldist2,Alltt2,Allrayp2,AlldistB,AllttB,AllraypB,&
     &lengthB,bitB,traveltimes,deltaf,synth_phase,synth_pick)


       use reloc_vars
       
       implicit none

       
	   real(8)      ALLdist2(nump0),ALLtt2(nump0), ALLrayp2(nump0)
	   integer      indexB,dirB,bitB, flagPASS
	   real(8)      AlldistB(nump0,maxB), AllttB(nump0,maxB), AllraypB(nump0,maxB)
	   real(8)      AlltempB(nump0)
	   integer      lengthB(maxB)
	   real(8)      checkNAN
	   integer      check99   
	   integer      mindistindex, KK, minraypindex
	   real(8)      mindistval  
       real(8)      FINALrayp, FINALdist, FINALtt,FINALraypB(maxB),FINALttB(maxB)
       integer      II, PP, WW, QQ
       real(8)      vel(*)
       integer      phase
       real(8)      traveltimes(*)
       real         epidist,deltaf(*)
       integer      synth_phase(synth_nq * synth_numstation)
       integer      synth_pick
       

     	   

       if ((phase.ne.1).and.(phase.ne.2)) THEN
    	   print*, 'Phase must be 1 (p) or 2 (S), now is:',phase
           stop
       end if
       

	   dirB = int((ALLdist2(2)-Alldist2(1))/abs(ALLdist2(2)-Alldist2(1)))
	   bitB = 1
	   indexB = 3

	   AlldistB = -9999
	   AllttB = -9999
	   AllraypB = -9999

	   AlldistB(1,1) = Alldist2(1)
	   AllttB(1,1) = Alltt2(1)
	   AllraypB(1,1) = Allrayp2(1)
	   AlldistB(2,1) = Alldist2(2)
	   AllttB(2,1) = Alltt2(2)
	   AllraypB(2,1) = Allrayp2(2)

	   DO II = 3,nump-1
	   
! 	     print*, Alldist2(II), alltt2(II), allrayp2(II), 1/vel(1)


		 IF (int((ALLdist2(II)-Alldist2(II-1))/abs(ALLdist2(II)-Alldist2(II-1))).ne.dirB) THEN

		   IF (dirB.eq.-1) THEN  !Reverse bitB to get increasing distances

			 AlldistB(1:indexB-3,bitB) = AlldistB(indexB-2:2:-1,bitB)
			 AllttB(1:indexB-3,bitB) = AllttB(indexB-2:2:-1,bitB)
			 AllraypB(1:indexB-3,bitB) = AllraypB(indexB-2:2:-1,bitB)

		   ELSE

			 AlldistB(1:indexB-3,bitB) = AlldistB(2:indexB-2,bitB)
			 AllttB(1:indexB-3,bitB) = AllttB(2:indexB-2,bitB)
			 AllraypB(1:indexB-3,bitB) = AllraypB(2:indexB-2,bitB)

		   END IF

		   lengthB(bitB) = indexB-3

		   checkNAN = 0.
		   check99 = 0
		   DO KK = 1,lengthB(bitB)
			 checkNAN = checkNAN + AlldistB(KK,bitB)
			 if (AlldistB(KK,bitB).eq.-9999) THEN
			   check99 = 1;
			 end if   
		   END DO

		   IF ((isnan(checkNAN)).or.(check99.eq.1)) THEN
		   ELSE     
			 IF (checkNAN.ne.0) THEN
			   bitB = bitB + 1
			 END IF         
		   END IF


		   dirB = int((ALLdist2(II)-Alldist2(II-1))/abs(ALLdist2(II)-Alldist2(II-1)))
		   indexB = 3

		   AlldistB(1,bitB) = Alldist2(II-1)
		   AllttB(1,bitB) = Alltt2(II-1)
		   AllraypB(1,bitB) = Allrayp2(II-1)
		   AlldistB(2,bitB) = Alldist2(II)
		   AllttB(2,bitB) = Alltt2(II)
		   AllraypB(2,bitB) = Allrayp2(II)
				   
		 ELSE

		   AlldistB(indexB,bitB) = Alldist2(II)
		   AllttB(indexB,bitB) = Alltt2(II)
		   AllraypB(indexB,bitB) = Allrayp2(II)

		   indexB = indexB + 1


		 END IF


	   END DO


	   IF (dirB.eq.-1) THEN  !Reverse bitB to get increasing distances

		 AlldistB(1:indexB-3,bitB) = AlldistB(indexB-2:2:-1,bitB)
		 AllttB(1:indexB-3,bitB) = AllttB(indexB-2:2:-1,bitB)
		 AllraypB(1:indexB-3,bitB) = AllraypB(indexB-2:2:-1,bitB)

	   END IF

	   lengthB(bitB) = indexB-3

	   IF (lengthB(bitB).eq.0) THEN   !Check if length is 0
		 bitB = bitB -1
	   END IF

	   IF (lengthB(bitB).eq.2) THEN   !Check if length is 2 and two are the same
		 IF ((AlldistB(1,bitB) + AlldistB(2,bitB)).eq.0) THEN
		   bitB = bitB -1
		 END IF
	   END IF
	   
! 	   print*,bitB,lengthB(bitB),'HAHA'

	   mindistindex = 1
	   mindistval = AlldistB(1,1)

	   !!  FIX short distances 0deg
	   IF (vel(2).gt.vel(1)) THEN     !Only do if gradient is positif at surface
		 DO II = 2,bitB
		   IF (AlldistB(1,II).lt.AlldistB(1,mindistindex)) THEN
			 mindistindex = II
			 mindistval = AlldistB(1,II)
		   END IF
		 END DO

		 !Distance   
		 AlltempB = AlldistB(:,mindistindex)
		 AlldistB(1,mindistindex) = 0.
		 AlldistB(2:lengthB(mindistindex)+1,mindistindex) = AlltempB(1:lengthB(mindistindex))
		 !TT
		 AlltempB = AllttB(:,mindistindex)
		 AllttB(1,mindistindex) = 0.
		 AllttB(2:lengthB(mindistindex)+1,mindistindex) = AlltempB(1:lengthB(mindistindex))
		 !rayp
		 AlltempB = AllraypB(:,mindistindex)
		 AllraypB(1,mindistindex) = 1/vel(1)     !pmax
		 AllraypB(2:lengthB(mindistindex)+1,mindistindex) = AlltempB(1:lengthB(mindistindex))   

		 lengthB(mindistindex) = lengthB(mindistindex) + 1

	   END IF
	   
	   
	   DO PP = 1,synth_pick
	   
	   FINALraypB = 1e9
	   FINALttB = 1e9	   
	   
! 		print*, 'C',synth_phase(PP),phase

	       IF (synth_phase(PP).eq.phase) THEN
	       

	       
           epidist = deltaf(PP)
	   
		   DO II = 1,bitB
		  
	!          print*,II,ALLdistB(1,II),ALLdistB(lengthB(II),II)
		  
			  IF ((epidist.ge.ALLdistB(1,II)).and.(epidist.le.ALLdistB(lengthB(II),II))) THEN
	!          print*,'JAJAJ'
				call interp_linear(1,lengthB(II),ALLdistB(1:lengthB(II),II),ALLraypB(1:lengthB(II),II), &
							  1,real(epidist,8),FINALraypB(II))
				call interp_linear(1,lengthB(II),ALLdistB(1:lengthB(II),II),ALLttB(1:lengthB(II),II),  &
							  1,real(epidist,8),FINALttB(II))
			  END IF
			
	!          print*,II,FINALraypB(II),FINALttB(II)
			END DO

!              if (PP.eq.79) then
!              	   DO QQ = 1,bitB
! 				 print*, lengthB(QQ)
! 				   do WW = 1,lengthB(QQ)
! 					write(100+PP,*) WW, ALLdistB(WW,QQ), ALLttB(WW,QQ)
! 				 end do		
! 				end do
!                 print*,'AHAHAHAHA',FINALttB(1:bitB), epidist
!              end if
              

			  IF (interpMETHOD.eq.1) THEN !minimize rayp
				  minraypINDEX =  minloc(FINALraypB(1:bitB),1)
				  FINALtt = FINALttB(minraypINDEX)
			  ELSE IF (interpMETHOD.eq.2) THEN !minimize tt
				  FINALtt = minval(FINALttB(1:bitB))
			  END IF
		  

		 
			 traveltimes(PP) = FINALtt
			 
			  write(*,*)   PP,deltaf(PP),synth_phase(PP),traveltimes(PP)
	    
	         END IF
       END DO
	   
! 	   print*,'D',synth_phase
	   
   
!        DO II = 1,bitB
!          print*,'bit:',II,lengthB(II)
!        END DO	   


END SUBROUTINE break_ttcurve_synth

SUBROUTINE read_file_picks

	use reloc_vars
	
      open (11,file=file_picks,status='old')
	
	  ipick=0
      do iq=1,nq0
         read (11,7,end=20) qid(iq),qlat(iq),qlon(iq),qdep(iq),&
     &           sc,sc_err,npick,etype(iq)
7        format (a10,1x,3f8.2,2f10.3,1x,i3,2x,a7)
         ipick1(iq)=ipick+1
         ipick2(iq)=ipick+npick
         do j=1,npick
            ipick=ipick+1
            read (11,*) istat,idph(ipick),pick(ipick),pickerr(ipick)
            if (istat.eq.12) then
               is=1
            else if (istat.eq.14) then
               is=2
            else if (istat.eq.15) then
               is=3
            else
               is=4
            end if
            idsta(ipick) = is
            qindex(ipick) = iq
            nstpick(is,1) = nstpick(is,1) + 1
         enddo
      enddo
      
20    close (11)

      nq=iq-1
      
!       print*, 'RFP:',nq
      
!       do iq = 1,ipick
!       print*,iq,idph(iq)
!       end do
      
END SUBROUTINE read_file_picks  

SUBROUTINE write_synth_data(vp,vs,Pdeltab8,Ptttab8,Pptab8,Sdeltab8,Stttab8,Sptab8)

	   use mt19937
	   use reloc_vars

	   implicit none

	   real        sta_lat(synth_numstation), sta_lon(synth_numstation)
	   real        epidist_synth(synth_nq * synth_numstation)
	   real(8)        Ptraveltimes(synth_nq * synth_numstation)
   	   real(8)        Straveltimes(synth_nq * synth_numstation)
	   integer     Psynth_phase(synth_nq * synth_numstation)
  	   integer     Ssynth_phase(synth_nq * synth_numstation)
	   integer     synth_pick, iq, kk,j
       real        lat_synth(synth_nq), lon_synth(synth_nq)

	   real(8)      PAlldistB(nump0,maxB), PAllttB(nump0,maxB), PAllraypB(nump0,maxB)
	   real(8)      SAlldistB(nump0,maxB), SAllttB(nump0,maxB), SAllraypB(nump0,maxB)  
	   real(8)      Pptab8(nump0), Ptttab8(nump0), Pdeltab8(nump0),Pztab8(nump0)
	   real(8)      Sptab8(nump0), Stttab8(nump0), Sdeltab8(nump0),Sztab8(nump0)
       integer      PlengthB(maxB), PbitB    
       integer      SlengthB(maxB), SbitB 
       
	   real(8)        vp(*),vs(*)	   

	   open (11,file='synthetic_data_input.txt',status='unknown')

	   sta_lat(1) =  -3.04
	   sta_lon(1) = -23.42
	   sta_lat(2) =  -3.65
	   sta_lon(2) = -17.48
	   sta_lat(3) =  26.08
	   sta_lon(3) =   3.66
	   sta_lat(4) =  -8.97
	   sta_lon(4) =  15.51

	   ! Compute random station locations if more than 4 apollo.
	   if (synth_numstation > 4) then
		  do kk = 5, synth_numstation - 4
   
			 sta_lat(4) =  genrand_real1() * 160 - 80   ! Latitude between -80 and 80
			 sta_lon(4) =  genrand_real1() * 160 - 80   ! Latitude between -80 and 80
   
		  end do
	   end if

	   do kk = 1,synth_nq

		 lat_synth(kk) = genrand_real1() * 160 - 80   ! Latitude between -80 and 80
		 lon_synth(kk) = genrand_real1() * 160 - 80   ! Longitude between -80 and 80

		 call SPH_DISTDP(lat_synth(kk),lon_synth(kk),sta_lat(1),sta_lon(1),epidist_synth((kk-1)*4 + 1))
		 call SPH_DISTDP(lat_synth(kk),lon_synth(kk),sta_lat(2),sta_lon(2),epidist_synth((kk-1)*4 + 2))
		 call SPH_DISTDP(lat_synth(kk),lon_synth(kk),sta_lat(3),sta_lon(3),epidist_synth((kk-1)*4 + 3))
		 call SPH_DISTDP(lat_synth(kk),lon_synth(kk),sta_lat(4),sta_lon(4),epidist_synth((kk-1)*4 + 4))
  
! 		Psynth_phase((kk-1)*4+1:(kk-1)*4+4) = 1   !P
! 		Ssynth_phase((kk-1)*4+1:(kk-1)*4+4) = 2   !S
		 
 
	   end do
	   
		Psynth_phase = 1   !P
		Ssynth_phase = 2   !S	   

	   synth_pick = synth_nq * synth_numstation
	   print*,'SP:',synth_pick

		   call break_ttcurve_synth(1,vp,Pdeltab8,Ptttab8,Pptab8,PAlldistB,&
		  &PAllttB,PAllraypB,PlengthB,PbitB,Ptraveltimes,epidist_synth,Psynth_phase,synth_pick)

		   call break_ttcurve_synth(2,vs,Sdeltab8,Stttab8,Sptab8,SAlldistB,&
		  &SAllttB,SAllraypB,SlengthB,SbitB,Straveltimes,epidist_synth,Ssynth_phase,synth_pick)  

	   ipick = 0
	   do iq= 1, synth_nq

		  write (11,7) iq,lat_synth(iq),lon_synth(iq),0.,&
	   &           0.,0.,synth_numstation*2,'SYN_NI'
7        format (i10.0,1x,3f8.2,2f10.3,1x,i3,2x,a7)
!          ipick1(iq)=ipick+1
!          ipick2(iq)=ipick+npick

		  do j=1,synth_numstation
			 if (j.eq.1) then
				is=12
			 else if (j.eq.2) then
				is=14
			 else if (j.eq.3) then
				is=15
			 else if (j.eq.4) then
				is=16
			 else
				is = j
			 end if
			 ipick=ipick+1

			 write(11,*) is,Psynth_phase(ipick),Ptraveltimes(ipick),genrand_real1()*4+1
! 			 ipick=ipick+1
			 write(11,*) is,Ssynth_phase(ipick),Straveltimes(ipick),genrand_real1()*4+1
!           
!             idsta(ipick) = is
!             qindex(ipick) = iq
!             nstpick(is,1) = nstpick(is,1) + 1
		  enddo
	   enddo

      close (11)

      
END SUBROUTINE write_synth_data  

