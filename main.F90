! =============================================================================
! This file is part of G2C3 version 1.0 
! G2C3 version 1.0 is released under the 3-Clause BSD license:
!
! Copyright (c) 2024, G2C3 Team (team leader: Animesh Kuley, akuley@iisc.ac.in)
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, 
!    this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, 
!    this list of conditions and the following disclaimer in the documentation 
!    and/or other materials provided with the distribution.
!
! 3. Neither the name of the GTC Team nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without 
!    specific prior written permission.
! ===============================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!       Global Gyrokinetic Code using Cylindriacl Coordinates (G2C3)         !
!                          Version 1.0                                       !
!                              2024                                          !
!                           G2C3 Team                                        !
!                                                                            !
! Developers:                                                                !
! Animesh Kuley (akuley@iisc.ac.in), , Jaya Kumar Alageshan,                 !
! Sarveshwar Sharma, Joydeep Das, Tajinder, Amal Biju, Vibhor Kumar Singh,   !
! Abhishek Tiwari, Shivam Prakash, Mudit                                     !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!====================================================================
program G2C3

use system_env
use global_parameters
use initial
use particle_array
use locator
use cylindricalRZ
use field_array
use Maps

implicit none

integer :: nnn,fid,i,istatus(MPI_STATUS_SIZE),ierror,j,ff,DataWriteStep,&
           k,l
real :: P, T, ss, NoiseAmplitude
real, allocatable, dimension(:) :: delN, PhiMyPP, PhiNextPP, Noise
real, allocatable, dimension(:,:,:) :: Phi3D
complex, allocatable, dimension(:,:,:) :: Spectrum
character*35 fname


SMOOTHFLAG = 1
PARALLELSMOOTHFLAG = 0
FILTERFLAG = 0
PROFILEFLAG = 1

FullFFLAG = 0

call G2C3_initial

call setup
write(*,*)'Setup done!'

call field_lines

call loadi_thermal

!!-----------------------------------------

! call Testing

!!--------------------------------------------

call random_number( zion(5,:) )
zion(5,:) = ( 0.5 - zion(5,:) ) * 0.1
zion(5,:) = 0.01 * zion(5,:)


!---------------------------------------------

! if( MyPE==0 )  then
!    write(fname,'("ZION",i2,".out")')(mype+10)
!    open(45+MyPE,file=fname,status='replace')
! 
!    do i=1,mi
!       write(45+MyPE,*)zion(:,i)
!    enddo
! 
!    close(45+MyPE)
! endif

!------------------------------------------------------

   allocate( PhiMyPP(mgridcore) )
   allocate( PhiNextPP(mgridcore) )
   allocate( Phi3D( NPsiBoozerGrid, NThetafBoozerGrid, NZetaBoozerGrid) )
   allocate( Spectrum( NPsiBoozerGrid , NThetafBoozerGrid , mtoroidal * CNZetaBoozerGrid ) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FEMphi = 0
FEMer = 0
FEMez = 0
FEMetor = 0

!!!!!!!!!!!!!!!!!!!


DataWriteStep = 100

do i=1,9000

   Iteration = i

   if( MyPE==0 )  call system("date |& tee -a Time.out")

   do irk=1,2

     call pushigc

     if( irk==2 )  call check_zone  ! ONLY for irk=2

     call shifti

     call FindGyroWeights
     call ScatterGyroDensity

     if( SMOOTHFLAG == 1 )  call SmoothGridField( FEMni )

     FEMni = FEMni * BCMask

     call Run_GK_Poisson_Solver

   enddo

   !--------------------------

   if( i==1 ) then
      call DumpInitData
   endif

   !--------------------------

  if( (i==1 .OR. modulo(i,DataWriteStep)==0) .AND. (modulo(MyPE,4)==0) ) then

    write(fname,'("PhiResult",i07,".out")')((mype+1)*10000+i)
    open(145+mype,file=fname,status='replace')

    do j=1,mgridcore
       write(145+mype,*)FEMphi(j)
    enddo

    close(145+mype)

    !---------------------------------
    write(fname,'("RhoResult",i07,".out")')((mype+1)*10000+i)
    open(145+mype,file=fname,status='replace')

    do j=1,mgridcore
       write(145+mype,*)FEMni(j)
    enddo

    close(145+mype)

    !---------------------------------

!    if( i>1000 ) then
!
!      write(fname,'("ER",i07,".out")')((mype+1)*10000+i)
!      open(145+mype,file=fname,status='replace')
!
!      do j=1,mgridcore
!         write(145+mype,*)FEMeR(j)
!      enddo
!    
!      close(145+mype)
!
!      !---------------------
!
!      write(fname,'("EZ",i07,".out")')((mype+1)*10000+i)
!      open(145+mype,file=fname,status='replace')
!
!      do j=1,mgridcore
!         write(145+mype,*)FEMeZ(j)
!      enddo
!    
!      close(145+mype)
!  
!      !---------------------
!
!      write(fname,'("EP",i07,".out")')((mype+1)*10000+i)
!      open(145+mype,file=fname,status='replace')
!
!      do j=1,mgridcore
!         write(145+mype,*)FEMeP(j)
!      enddo
!    
!      close(145+mype)
!
!      !---------------------
!
!      write(fname,'("ET",i07,".out")')((mype+1)*10000+i)
!      open(145+mype,file=fname,status='replace')
!
!      do j=1,mgridcore
!         write(145+mype,*)FEMetor(j)
!      enddo
!    
!      close(145+mype)
!
!
!    endif

!    if( i>2900 ) then
!
!      write(fname,'("zzion",i07,".out")')((mype+1)*10000+i)
!      open(145+mype,file=fname,status='replace')
!
!      do j=1,mi
!         write(145+mype,*)zion(:,j)
!      enddo
!    
!      close(145+mype)
!
!    endif

    if( i>1000 ) then

      call Field2BoozerGridField( FEMPhi , Phi3D )

      write(fname,'("Phi3D",i07,".out")')((mype+1)*10000+i)
      open(145+mype,file=fname,status='replace')

      write(145+mype,*)NPsiBoozerGrid
      write(145+mype,*)NThetafBoozerGrid
      write(145+mype,*)NZetaBoozerGrid

      do j=1,NPsiBoozerGrid
        do k=1,NThetafBoozerGrid
          do l=1,NZetaBoozerGrid
            write(145+mype,*)Phi3D(j,k,l)
          enddo
        enddo
      enddo
       
      close(145+mype)

!      !-----------------------------------
!
!      Phi3D = 0.0
!
!      call Field2BoozerGridField( FEMeR , Phi3D )
!
!      write(fname,'("ER3D",i07,".out")')((mype+1)*10000+i)
!      open(145+mype,file=fname,status='replace')
!
!      write(145+mype,*)NPsiBoozerGrid
!      write(145+mype,*)NThetafBoozerGrid
!      write(145+mype,*)NZetaBoozerGrid
!
!      do j=1,NPsiBoozerGrid
!        do k=1,NThetafBoozerGrid
!          do l=1,NZetaBoozerGrid
!            write(145+mype,*)Phi3D(j,k,l)
!          enddo
!        enddo
!      enddo
!
!      close(145+mype)
!
!      !-----------------------------------
!
!      Phi3D = 0.0
!
!      call Field2BoozerGridField( FEMeZ , Phi3D )
!
!      write(fname,'("EZ3D",i07,".out")')((mype+1)*10000+i)
!      open(145+mype,file=fname,status='replace')
!
!      write(145+mype,*)NPsiBoozerGrid
!      write(145+mype,*)NThetafBoozerGrid
!      write(145+mype,*)NZetaBoozerGrid
!
!      do j=1,NPsiBoozerGrid
!        do k=1,NThetafBoozerGrid
!          do l=1,NZetaBoozerGrid
!            write(145+mype,*)Phi3D(j,k,l)
!          enddo
!        enddo
!      enddo
!
!      close(145+mype)
!
!      !-----------------------------------
!
!      Phi3D = 0.0
!
!      call Field2BoozerGridField( FEMeTor , Phi3D )
!
!      write(fname,'("ET3D",i07,".out")')((mype+1)*10000+i)
!      open(145+mype,file=fname,status='replace')
!
!      write(145+mype,*)NPsiBoozerGrid
!      write(145+mype,*)NThetafBoozerGrid
!      write(145+mype,*)NZetaBoozerGrid
!
!      do j=1,NPsiBoozerGrid
!        do k=1,NThetafBoozerGrid
!          do l=1,NZetaBoozerGrid
!            write(145+mype,*)Phi3D(j,k,l)
!          enddo
!        enddo
!      enddo
!
!      close(145+mype)
!
!      !-----------------------------------

    endif

  endif

enddo

!---------------------------------

if( MyPE==0 ) then
   write(fname,'("zzion",i2,".out")')(mype+10)
   open(45+mype,file=fname,status='replace')

   do j=1,mi
      write(45+mype,*)zion(:,j)
   enddo

   close(45+mype)
endif

!---------------------------------

if( MyPE==0 )  then
   write(G2C3out,*)"Reached end of main.F90...!"
endif

flush(G2C3out)

call G2C3_final

end program G2C3
