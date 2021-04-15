!---------------------------------------------------------------------!
! Created by Madu Manathunga on 06/29/2020                            !
!                                                                     ! 
! Previous contributors: Yipu Miao, Xio He, Alessandro Genoni,        !
!                         Ken Ayers & Ed Brothers                     !
!                                                                     !
! Copyright (C) 2020-2021 Merz lab                                    !
! Copyright (C) 2020-2021 Götz lab                                    !
!                                                                     !
! This Source Code Form is subject to the terms of the Mozilla Public !
! License, v. 2.0. If a copy of the MPL was not distributed with this !
! file, You can obtain one at http://mozilla.org/MPL/2.0/.            !
!_____________________________________________________________________!

!---------------------------------------------------------------------!
! This module contains subroutines and data structures related to     ! 
! scf gradient calculation.                                           !
!---------------------------------------------------------------------!

#include "util.fh"

#ifdef OSHELL
module quick_oshell_gradient_module
#else
module quick_cshell_gradient_module
#endif

  implicit double precision(a-h,o-z)
  private

#ifdef OSHELL
   public  :: oshell_gradient, uscf_gradient, get_oshell_oneen_grad
#else
  public  :: cshell_gradient, get_cshell_oneen_grad, scf_gradient, get_nuclear_repulsion_grad, get_ijbas_derivative
  public :: allocate_quick_gradient, deallocate_quick_gradient, tmp_grad, tmp_ptchg_grad
  
  double precision, allocatable, dimension(:) :: tmp_grad
  double precision, allocatable, dimension(:) :: tmp_ptchg_grad
#endif

contains

#ifndef OSHELL
  subroutine allocate_quick_gradient()

    use quick_molspec_module, only : quick_molspec, natom
    implicit none

    if(.not. allocated(tmp_grad)) allocate(tmp_grad(3*natom))
    tmp_grad = 0.0d0

    if(quick_molspec%nextatom.gt.0) then
      if(.not. allocated(tmp_ptchg_grad)) allocate(tmp_ptchg_grad(3*quick_molspec%nextatom))    
      tmp_ptchg_grad = 0.0d0
    endif

  end subroutine allocate_quick_gradient

  subroutine deallocate_quick_gradient()

    use quick_molspec_module, only : quick_molspec
    implicit none

    if(allocated(tmp_grad)) deallocate(tmp_grad) 
    
    if(quick_molspec%nextatom.gt.0) then
      if(allocated(tmp_ptchg_grad)) deallocate(tmp_ptchg_grad)
    endif
  
  end subroutine deallocate_quick_gradient
#endif


#ifdef OSHELL 
  subroutine oshell_gradient(ierr)
#else
  subroutine cshell_gradient(ierr)
#endif
  
  !------------------------------------------------------------------
  ! This subroutine carries out a gradient calculation 
  !------------------------------------------------------------------
  
     use allmod
     implicit double precision(a-h,o-z)
  
     integer, intent(inout) :: ierr
     character(len=1) cartsym(3)
  
#ifdef MPIV
     include "mpif.h"
#endif
  
  !  Curently, only analytical gradients are available. This should be changed later.
     quick_method%analgrad=.true.
  
     quick_method%integralCutoff=1.0d0/(10.0d0**6.0d0)
     quick_method%Primlimit=1.0d0/(10.0d0**6.0d0)
  
  !  Set array elements required for printing the gradients
     cartsym(1) = 'X'
     cartsym(2) = 'Y'
     cartsym(3) = 'Z'
  
  !  Set the value of gradient vector to zero
     do j=1,natom
        do k=1,3
           quick_qm_struct%gradient((j-1)*3+K)=0d0
        enddo
     enddo
  
     call getEnergy(.false.,ierr)
  
     if (quick_method%analgrad) then
#ifdef OSHELL
        call uscf_gradient
#else
        call scf_gradient
#endif
     endif
  
#if defined CUDA || defined CUDA_MPIV
     if (quick_method%bCUDA) then
        call gpu_cleanup()
     endif
#endif
  
#ifdef MPIV
     if(master) then
#endif
  
     call PrtAct(ioutfile,"Begin Gradient Calculation")
     write (ioutfile,'(" ANALYTICAL GRADIENT: ")')
     write (ioutfile,'(40("-"))')
     write (ioutfile,'(" COORDINATE",4x,"XYZ",12x,"GRADIENT")')
     write (ioutfile,'(40("-"))')
     do Iatm=1,natom
        do Imomentum=1,3
           write (ioutfile,'(I5,A1,3x,F14.10,3x,F14.10)')Iatm,cartsym(imomentum), &
           xyz(Imomentum,Iatm)*0.529177249d0,quick_qm_struct%gradient((Iatm-1)*3+Imomentum)
        enddo
     enddo
  
     write(ioutfile,'(40("-"))')
     
     if(quick_method%extCharges) then
        write (ioutfile,'(/," POINT CHARGE GRADIENT: ")')
        write (ioutfile,'(40("-"))')
        write (ioutfile,'(" COORDINATE",4x,"XYZ",12x,"GRADIENT")')
        write (ioutfile,'(40("-"))')
        do Iatm=1,quick_molspec%nextatom
           do Imomentum=1,3
              write (ioutfile,'(I5,A1,3x,F14.10,3x,F14.10)')Iatm,cartsym(imomentum), &
              quick_molspec%extxyz(Imomentum,Iatm)*0.529177249d0,quick_qm_struct%ptchg_gradient((Iatm-1)*3+Imomentum)
           enddo
        enddo
        write(ioutfile,'(40("-"))')
     endif   
  
     call PrtAct(ioutfile,"End Gradient Calculation")
  
#ifdef MPIV
     endif
#endif
  
     return

#ifdef OSHELL
  end subroutine oshell_gradient
#else
  end subroutine cshell_gradient
#endif
  
  
#ifdef OSHELL
  subroutine uscf_gradient
#else
  subroutine scf_gradient
#endif
     use allmod
#ifdef OSHELL
     use quick_cshell_gradient_module
#endif
     implicit double precision(a-h,o-z)
  
     integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
     common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
#ifdef MPIV
     include "mpif.h"
#endif
  
  !---------------------------------------------------------------------
  !  The purpose of this subroutine is to calculate the gradient of
  !  the total energy with respect to nuclear displacement.  The results
  !  of this are stored in Gradient, which is organized by atom and then
  !  by direction of displacement, i.e. element 1 is the gradient of the
  !  x diplacement of atom 1, element 5 is the y displacement of atom 2.
  !  Not that this is the RHF version of the code.  It is simplest of
  !  the gradient codes in this program.
  !  The gradient at this level of theory is the sum of five terms.
  !---------------------------------------------------------------------
  
  !  Start the timer for gradient calculation
     call cpu_time(timer_begin%TGrad)
  
#ifdef MPIV
     call allocate_quick_gradient()
#endif
  
  !  Set the values of gradient arry to zero 
     quick_qm_struct%gradient       = 0.0d0
     if (quick_molspec%nextatom .gt. 0) quick_qm_struct%ptchg_gradient = 0.0d0
  
  !---------------------------------------------------------------------
  !  1) The derivative of the nuclear repulsion.
  !---------------------------------------------------------------------
  
     call cpu_time(timer_begin%TNucGrad)
  
     call get_nuclear_repulsion_grad
  
     call cpu_time(timer_end%TNucGrad)
     timer_cumer%TNucGrad = timer_cumer%TNucGrad + timer_end%TNucGrad-timer_begin%TNucGrad
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!Madu!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPIV
  if(master) then
#endif
  
!#ifdef DEBUG
!     if (quick_method%debug) then
          write (*,'(/," DEBUG STEP 1 :  NUCLEAR REPULSION GRADIENT: ")')
          do Iatm=1,natom
              do Imomentum=1,3
                  write (*,'(I5,7x,F20.10)')Iatm, &
                  quick_qm_struct%gradient((Iatm-1)*3+Imomentum)
              enddo
          enddo
!     endif
!#endif
  
#ifdef MPIV
  endif
#endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!Madu!!!!!!!!!!!!!!!!!!!!!!!
  
  !---------------------------------------------------------------------
  !  2) One electron gradients
  !---------------------------------------------------------------------
  ! Note that we will call this subroutine asynchronously with ERI
  ! gradient kernel call (see gpu_get2e.cu) in CUDA and CUDA_MPI versions
  
!#if !defined CUDA && !defined CUDA_MPIV
#ifdef OSHELL
    call get_oshell_oneen_grad
#else
    call get_cshell_oneen_grad
#endif
!#endif
  !---------------------------------------------------------------------
  !  3) The derivative of the electron repulsion term
  !---------------------------------------------------------------------
#ifdef MPIV
     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
#endif
  
     call cpu_time(timer_begin%T2eGrad)
  
#ifdef OSHELL
     call get_oshell_eri_grad
#else
     call get_cshell_eri_grad
#endif
  
     call cpu_time(timer_end%T2eGrad)
     timer_cumer%T2eGrad = timer_cumer%T2eGrad + timer_end%T2eGrad-timer_begin%T2eGrad
  !---------------------------------------------------------------------
  !  4) If DFT, calculate the derivative of exchahnge correlation  term
  !---------------------------------------------------------------------
#ifdef MPIV
     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
#endif
  
     if (quick_method%DFT) then
        call cpu_time(timer_begin%TExGrad)
  
#ifdef OSHELL
        call get_oshell_xc_grad
#else
        call get_cshell_xc_grad
#endif
  
#ifdef CUDA_MPIV
        call mgpu_get_xcrb_time(timer_cumer%TDFTrb, timer_cumer%TDFTpg)
#endif
  
        call cpu_time(timer_end%TExGrad)
        timer_cumer%TExGrad = timer_cumer%TExGrad + timer_end%TExGrad-timer_begin%TExGrad
  
     endif
  
  
#ifdef MPIV
  
     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
  
     call cpu_time(timer_begin%TGradred) 
  
  ! sum up all gradient contributions
     call MPI_REDUCE(quick_qm_struct%gradient, tmp_grad, 3*natom, mpi_double_precision, MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
     if(quick_molspec%nextatom.gt.0) call MPI_REDUCE(quick_qm_struct%ptchg_gradient, tmp_ptchg_grad, 3*quick_molspec%nextatom,& 
                                     mpi_double_precision, MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
     if(master) then
       quick_qm_struct%gradient(:) = tmp_grad(:)
       if(quick_molspec%nextatom.gt.0) quick_qm_struct%ptchg_gradient(:) = tmp_ptchg_grad(:)
     endif
  
     call cpu_time(timer_end%TGradred)
  
     timer_cumer%TGradred = timer_cumer%TGradred + timer_end%TGradred-timer_begin%TGradred
  
#endif
  
  !!!!!!!!!!!!!!!!!!!!!!!Madu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPIV
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
  if(master) then
#endif
  
!#ifdef DEBUG
!    if (quick_method%debug) then
          write (*,'(/," DEBUG STEP : TOTAL GRADIENT: ")')
          do Iatm=1,natom*3
                  write (*,'(I5,7x,F20.10)')Iatm,quick_qm_struct%gradient(Iatm)
          enddo
!    endif
!#endif
  
#ifdef MPIV
  endif
#endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!Madu!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!!!!!Madu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPIV
  if(master) then
#endif
  
#ifdef DEBUG
    if (quick_method%debug) then
      if(quick_method%extCharges) then
          write (iOutFile,'(/," DEBUG: POINT CHARGE GRADIENT: ")')
          do Iatm=1,quick_molspec%nextatom*3
                  write (iOutFile,'(I5,7x,F20.10)')Iatm,quick_qm_struct%ptchg_gradient(Iatm)
          enddo
      endif
    endif
  
#endif
  
#ifdef MPIV
  endif
#endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!Madu!!!!!!!!!!!!!!!!!!!!!!!
  
#ifdef MPIV
     call deallocate_quick_gradient()
#endif
  
  !  Stop the timer and add up the total gradient times
     call cpu_time(timer_end%TGrad)
     timer_cumer%TGrad=timer_cumer%TGrad+timer_end%TGrad-timer_begin%TGrad
  
     return
  
#ifdef OSHELL
  end subroutine uscf_gradient
#else
  end subroutine scf_gradient
#endif
  
#ifndef OSHELL
  subroutine get_nuclear_repulsion_grad
  
     use allmod
     implicit double precision(a-h,o-z)
  
     double precision, external :: rootSquare
     integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
     common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  
  !  This subroutine calculates the nuclear repulsion gradients. 
  
  !  Quick derivation:
  !  Vnn = (Sum over A) (Sum over B>A) ZA ZB / RAB
  !  where A and B are atoms, Z are charges, and RAB is the interatomic
  !  seperation.  If we take the derivative, all terms not involving
  !  A fall out. Thus:
  !  Vnn/dXA = ZA (Sum over B) d/dXA (ZB /RAB)
  !  Vnn/dXA = ZA (Sum over B) ZB d/dXA (RAB^-1)
  !  Vnn/dXA = ZA (Sum over B) ZB d/dXA(((XA-XB)^2+(YA-YB)^2+(ZA-ZB)^2)^-.5)
  !  Vnn/dXA = ZA (Sum over B) ZB*-.5*((XA-XB)^2+(YA-YB)^2+(ZA-ZB)^2)^-1.5
  !  *2*(XA-XB)^1
  !  Vnn/dXA = ZA (Sum over B) ZB*-((XA-XB)^2+(YA-YB)^2+(ZA-ZB)^2)^-1.5
  !  *(XA-XB)^1
  !  Vnn/dXA = ZA (Sum over B) ZB*((XA-XB)^2+(YA-YB)^2+(ZA-ZB)^2)^-1.5
  !  *(XB-XA)
  !  Vnn/dXA = ZA (Sum over B) ZB*(XB-XA) RAB^-3
  !
  !  Thus this term is trivial, and can calculate it here. Note also
  !  that that atom A can never equal atom B, and A-B part of the derivative
  !  for A is the negative of the BA derivative for atom B.
  
     do Iatm = 1,(natom+quick_molspec%nextatom)*3
        do Jatm = Iatm+1,(natom+quick_molspec%nextatom)
           if(master) then
              if(Iatm<=natom .and. Jatm<=natom)then  
  !  Nuclear-nuclear repulsion grdients, update nuclear gradient vector
                 RIJ = rootSquare(xyz(1:3,Iatm), xyz(1:3,Jatm), 3)
                 ZAZBdivRIJ3 = quick_molspec%chg(Iatm)*quick_molspec%chg(Jatm)*(RIJ**(-3.0d0))
                 XBminXA = xyz(1,Jatm)-xyz(1,Iatm)
                 YBminYA = xyz(2,Jatm)-xyz(2,Iatm)
                 ZBminZA = xyz(3,Jatm)-xyz(3,Iatm)
                 ISTART = (Iatm-1)*3
                 JSTART = (Jatm-1)*3
                 quick_qm_struct%gradient(ISTART+1) = quick_qm_struct%gradient(ISTART+1)+XBminXA*ZAZBdivRIJ3
                 quick_qm_struct%gradient(ISTART+2) = quick_qm_struct%gradient(ISTART+2)+YBminYA*ZAZBdivRIJ3
                 quick_qm_struct%gradient(ISTART+3) = quick_qm_struct%gradient(ISTART+3)+ZBminZA*ZAZBdivRIJ3
                 quick_qm_struct%gradient(JSTART+1) = quick_qm_struct%gradient(JSTART+1)-XBminXA*ZAZBdivRIJ3
                 quick_qm_struct%gradient(JSTART+2) = quick_qm_struct%gradient(JSTART+2)-YBminYA*ZAZBdivRIJ3
                 quick_qm_struct%gradient(JSTART+3) = quick_qm_struct%gradient(JSTART+3)-ZBminZA*ZAZBdivRIJ3
  
              elseif(Iatm<=natom .and. Jatm>natom)then 
  
  !  Nuclear-point charge repulsion grdients, update nuclear gradient vector               
  
                 RIJ = rootSquare(xyz(1:3,Iatm), quick_molspec%extxyz(1:3,Jatm-natom), 3)
                 ZAZBdivRIJ3 = quick_molspec%chg(Iatm)*quick_molspec%extchg(Jatm-natom)*(RIJ**(-3.0d0))
                 XBminXA = quick_molspec%extxyz(1,Jatm-natom)-xyz(1,Iatm)
                 YBminYA = quick_molspec%extxyz(2,Jatm-natom)-xyz(2,Iatm)
                 ZBminZA = quick_molspec%extxyz(3,Jatm-natom)-xyz(3,Iatm) 
                 ISTART = (Iatm-1)*3
                 quick_qm_struct%gradient(ISTART+1) = quick_qm_struct%gradient(ISTART+1)+XBminXA*ZAZBdivRIJ3
                 quick_qm_struct%gradient(ISTART+2) = quick_qm_struct%gradient(ISTART+2)+YBminYA*ZAZBdivRIJ3
                 quick_qm_struct%gradient(ISTART+3) = quick_qm_struct%gradient(ISTART+3)+ZBminZA*ZAZBdivRIJ3
             
  !  Nuclear-point charge repulsion grdients, update point charge gradient vector
  
                 JSTART = (Jatm-natom-1)*3
                 quick_qm_struct%ptchg_gradient(JSTART+1) = quick_qm_struct%ptchg_gradient(JSTART+1)-XBminXA*ZAZBdivRIJ3
                 quick_qm_struct%ptchg_gradient(JSTART+2) = quick_qm_struct%ptchg_gradient(JSTART+2)-YBminYA*ZAZBdivRIJ3
                 quick_qm_struct%ptchg_gradient(JSTART+3) = quick_qm_struct%ptchg_gradient(JSTART+3)-ZBminZA*ZAZBdivRIJ3
              endif
  
           endif
        enddo
     enddo
  
     return
  
  end subroutine get_nuclear_repulsion_grad
#endif  

 
#ifdef OSHELL 
  subroutine get_oshell_oneen_grad
#else
  subroutine get_cshell_oneen_grad
#endif  

    use allmod
    use quick_oei_module, only:attrashellopt
    implicit none
    integer :: Iatm, Imomentum, IIsh, JJsh, i, j, nshell_mpi
  
#ifdef MPIV
     include "mpif.h"
#endif
  
  !---------------------------------------------------------------------
  !  1) The derivative of the kinetic term
  !---------------------------------------------------------------------
  
     call cpu_time(timer_begin%T1eGrad)
  
     call cpu_time(timer_begin%T1eTGrad)
  
#ifdef OSHELL
     call get_oshell_kinetic_grad
#else
     call get_cshell_kinetic_grad
#endif
  
     call cpu_time(timer_end%T1eTGrad)
     timer_cumer%T1eTGrad=timer_cumer%T1eTGrad+timer_end%T1eTGrad-timer_begin%T1eTGrad
  
#ifdef MPIV
  if(master) then
#endif
  
!#ifdef DEBUG
!    if (quick_method%debug) then
          write (*,'(/," DEBUG STEP 2 :  KINETIC GRADIENT ADDED: ")')
          do Iatm=1,natom
              do Imomentum=1,3
                  write (*,'(I5,7x,F20.10)')Iatm, &
                  quick_qm_struct%gradient((Iatm-1)*3+Imomentum)
              enddo
          enddo
!    endif
!#endif
  
#ifdef MPIV
  endif
#endif
  
  !---------------------------------------------------------------------
  !  2) The derivative of the 1 electron nuclear attraction term ij times
  !     the density matrix element ij.
  !---------------------------------------------------------------------
  
     call cpu_time(timer_begin%T1eVGrad)
  
#ifdef MPIV
     if (bMPI) then
        nshell_mpi = mpi_jshelln(mpirank)
     else
        nshell_mpi = jshell
     endif
  
     do i=1,nshell_mpi
        if (bMPI) then
           IIsh = mpi_jshell(mpirank,i)
        else
           IIsh = i
        endif
#else
     do IIsh=1,jshell
#endif
        do JJsh=IIsh,jshell
           call attrashellopt(IIsh,JJsh)
        enddo
     enddo
  
     call cpu_time(timer_end%T1eVGrad)
     timer_cumer%T1eVGrad=timer_cumer%T1eVGrad+timer_end%T1eVGrad-timer_begin%T1eVGrad
  
#ifdef MPIV
     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
#endif
  
#ifdef MPIV
  if(master) then
#endif
  
!#ifdef DEBUG
!    if (quick_method%debug) then
          write (*,'(/," DEBUG STEP 3 :  NUC-EN ATTRACTION GRADIENT ADDED:")')
          do Iatm=1,natom
              do Imomentum=1,3
                  write (*,'(I5,7x,F20.10)')Iatm, &
                  quick_qm_struct%gradient((Iatm-1)*3+Imomentum)
              enddo
          enddo
!    endif
!#endif
  
#ifdef MPIV
  endif
#endif
  
     call cpu_time(timer_end%T1eGrad)
     timer_cumer%T1eGrad = timer_cumer%T1eGrad + timer_end%T1eGrad-timer_begin%T1eGrad
  
     return
 
#ifdef OSHELL 
  end subroutine get_oshell_oneen_grad
#else
  end subroutine get_cshell_oneen_grad
#endif 
  
  
 
#ifdef OSHELL
  subroutine get_oshell_kinetic_grad
#else
  subroutine get_cshell_kinetic_grad
#endif 

     use allmod
#ifdef OSHELL
     use quick_cshell_gradient_module, only:get_ijbas_derivative
#endif

     implicit double precision(a-h,o-z)
  
     integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
     common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
     logical :: ijcon
#ifdef MPIV
     include "mpif.h"
#endif
  
  !  1)  The negative of the energy weighted density matrix element i j
  !   with the derivative of the ij overlap.
  !
  !  2)  The derivative of the 1 electron kinetic energy term ij times
  !  the density matrix element ij.
  !
  !  These terms are grouped together since we loop over the same terms.
  !  Also note that these are the 2-center terms.
  !
  !  The energy weighted denisty matrix is:
  !  Q(i,j) =2*(Sum over alpha electrons a)  E(a) C(I,a) C(J,a)
  !  Where C is the alpha or beta molecular orbital coefficients, and
  !  E is the alpha or beta molecular orbital energies.
  !  We'll store this in HOLD as we don't really need it (except for hessian
  !  calculations later).
  
  !write(*,*) "get_nuclear_repulsion_grad: Calculating HOLD array"
  
#ifdef MPIV
     if (master) then
#endif
     do I=1,nbasis
        do J=1,nbasis
           HOLDJI = 0.d0
#ifdef OSHELL
         do k=1,quick_molspec%nelec
            HOLDJI =HOLDJI+(quick_qm_struct%E(K)*quick_qm_struct%co(J,K)*quick_qm_struct%co(I,K))
         enddo

         do k=1,quick_molspec%nelecb
            HOLDJI =HOLDJI+(quick_qm_struct%Eb(K)*quick_qm_struct%cob(J,K)*quick_qm_struct%cob(I,K))
         enddo
         quick_scratch%hold(J,I)=HOLDJI

#else
           do K=1,quick_molspec%nelec/2
           HOLDJI = HOLDJI + (quick_qm_struct%E(K)*quick_qm_struct%co(J,K)*quick_qm_struct%co(I,K))
           enddo
           quick_scratch%hold(J,I) = 2.d0*HOLDJI
#endif
        enddo
     enddo
#ifdef MPIV
     endif
     call MPI_BCAST(quick_scratch%hold,nbasis*nbasis,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
#endif
  
     if (quick_method%debug) then
        write(ioutfile,'(/"THE ENERGY WEIGHTED DENSITY MATRIX")')
        do I=1,nbasis
           do J=1,nbasis
              write (ioutfile,'("W[",I4,",",I4,"]=",F18.10)') &
              J,I,quick_scratch%hold(J,I)
           enddo
        enddo
     endif
  
  !  The contribution to the derivative of energy with respect to nuclear
  !  postion for this term is: -(Sum over i,j) Q(i,j) dS(ij)/dXA
  !  Now Q is symmetric, and dS(ij)/dXA = dS(ji)/dXA.  Furthermore, if
  !  i and j are on the same center, the term is zero. Thus we need to find
  !  the i j pairs for i and j not on the same atom.
  !  Also:  The derivative of a cartesian gtf is:
  !
  !  d/dXA ((x-XA)^i (y-YA)^j (z-ZA)^k e^(-ar^2))
  !  = 2a((x-XA)^(i+1) (y-YA)^j (z-ZA)^k e^(-ar^2))
  !  - i ((x-XA)^(i-1) (y-YA)^j (z-ZA)^k e^(-ar^2))
  
  !  Note that the negative on the final term comes from the form of (x-XA).
#ifdef MPIV
     if (bMPI) then
        nbasis_mpi = mpi_nbasisn(mpirank)
     else
        nbasis_mpi = nbasis
     endif
  
     do i=1,nbasis_mpi
        if (bMPI) then
           Ibas = mpi_nbasis(mpirank,i)
        else
           Ibas = i
     endif
#else
     do Ibas=1,nbasis
#endif
        ISTART = (quick_basis%ncenter(Ibas)-1) *3
           do Jbas=quick_basis%last_basis_function(quick_basis%ncenter(IBAS))+1,nbasis
  
              JSTART = (quick_basis%ncenter(Jbas)-1) *3
#ifdef OSHELL
            DENSEJI = quick_qm_struct%dense(Jbas,Ibas)+quick_qm_struct%denseb(Jbas,Ibas)
#else
              DENSEJI = quick_qm_struct%dense(Jbas,Ibas)
#endif
  
  !  We have selected our two basis functions, now loop over angular momentum.
              do Imomentum=1,3
  
  !  do the Ibas derivatives first. In order to prevent code duplication,
  !  this has been implemented in a seperate subroutine. 
                 ijcon = .true. 
                 call get_ijbas_derivative(Imomentum, Ibas, Jbas, Ibas, ISTART, ijcon, DENSEJI) 
  
  !  do the Jbas derivatives.
                 ijcon = .false.
                 call get_ijbas_derivative(Imomentum, Ibas, Jbas, Jbas, JSTART, ijcon, DENSEJI)
  
              enddo
           enddo
        enddo
  
#ifdef MPIV
     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
#endif
  
     return
  
#ifdef OSHELL
  end subroutine get_oshell_kinetic_grad
#else
  end subroutine get_cshell_kinetic_grad
#endif
  
  
#ifdef OSHELL
  subroutine get_oshell_eri_grad
#else
  subroutine get_cshell_eri_grad
#endif
  
     use allmod
#ifdef OSHELL
   use quick_cutoff_module, only: oshell_density_cutoff, oshell_dnscreen
   use quick_oshell_eri_grad_module, only: oshell_eri_grad
#else
   use quick_cutoff_module, only:cshell_density_cutoff, cshell_dnscreen
   use quick_cshell_eri_grad_module, only: cshell_eri_grad
#endif
     implicit double precision(a-h,o-z)
  
     integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
     common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
    
#ifdef MPIV
     include "mpif.h"
#endif
  
  !  This subroutine calculates the derivative of 4center 2e- 
  !  terms with respect to X times the coefficient found in the energy.
  !  (i.e. the multiplicative constants from the density matrix that 
  !  arise as these are both the exchange and correlation integrals.

#ifdef OSHELL
      call oshell_density_cutoff
#else
      call cshell_density_cutoff
#endif
  
#if defined CUDA || defined CUDA_MPIV
     if (quick_method%bCUDA) then
  
        call gpu_upload_density_matrix(quick_qm_struct%dense)
        call gpu_upload_cutoff(cutmatrix, quick_method%integralCutoff,quick_method%primLimit,quick_method%DMCutoff)
        call gpu_upload_grad(quick_method%gradCutoff)

#ifdef OSHELL
        call gpu_upload_calculated_beta(quick_qm_struct%ob,quick_qm_struct%denseb)
        call gpu_get_oshell_eri_grad(quick_qm_struct%gradient)
#else
        call gpu_get_cshell_eri_grad(quick_qm_struct%gradient)
#endif
  
     else
#endif
  
#if defined MPIV && !defined CUDA_MPIV 
  
     if (bMPI) then
        nshell_mpi = mpi_jshelln(mpirank)
     else
        nshell_mpi = jshell
     endif
  
     do i=1,nshell_mpi
        if (bMPI) then
           II = mpi_jshell(mpirank,i)
        else
           II = i
     endif
#else
        do II=1,jshell
#endif
           do JJ=II,jshell
           Testtmp=Ycutoff(II,JJ)
              do KK=II,jshell
                 do LL=KK,jshell
                    if(quick_basis%katom(II).eq.quick_basis%katom(JJ).and.quick_basis%katom(II).eq. &
                    quick_basis%katom(KK).and.quick_basis%katom(II).eq.quick_basis%katom(LL))then
                       continue
                    else
                       testCutoff = TESTtmp*Ycutoff(KK,LL)
                       if(testCutoff.gt.quick_method%gradCutoff)then
                          DNmax=max(4.0d0*cutmatrix(II,JJ),4.0d0*cutmatrix(KK,LL), &
                          cutmatrix(II,LL),cutmatrix(II,KK),cutmatrix(JJ,KK),cutmatrix(JJ,LL))
                          cutoffTest=testCutoff*DNmax
                          if(cutoffTest.gt.quick_method%gradCutoff)then
#ifdef OSHELL
                             call oshell_eri_grad
#else
                             call cshell_eri_grad
#endif
                          endif
                       endif
                    endif
                 enddo
              enddo
           enddo
        enddo
#if defined CUDA || defined CUDA_MPIV
     endif
#endif

!#ifdef DEBUG
!    if (quick_method%debug) then
          write (*,'(/," DEBUG STEP 4 :  ERI GRADIENT ADDED:")')
          do Iatm=1,natom
              do Imomentum=1,3
                  write (*,'(I5,7x,F20.10)')Iatm, &
                  quick_qm_struct%gradient((Iatm-1)*3+Imomentum)
              enddo
          enddo
!    endif
!#endif
  
#ifdef MPIV
     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
#endif
  
     return
  
#ifdef OSHELL
  end subroutine get_oshell_eri_grad
#else
  end subroutine get_cshell_eri_grad
#endif
  
  
#ifdef OSHELL
  subroutine get_oshell_xc_grad
#else
  subroutine get_cshell_xc_grad
#endif
  
  !-------------------------------------------------------------------------
  !  This subroutine will calculate:
  !  1) The derivative of the exchange/correlation functional energy
  !  with respect to nuclear displacement.
  !  2) The derivative of the weight of the quadrature points with respect
  !  to nuclear displacement.
  !
  !  These two terms arise because of the quadrature used to calculate the
  !  XC terms.
  !  Exc = (Sum over grid points) W(g) f(g)
  !  dExc/dXA = (Sum over grid points) dW(g)/dXA f(g) + W(g) df(g)/dXA
  !
  !  For the W(g) df(g)/dXA term, the derivation was done by Ed Brothers and
  !  is a varient of the method found in the Johnson-Gill-Pople paper.  It can
  !  be found in Ed's thesis, assuming he ever writes it.
  !
  !  One of the actuals element is:
  !  dExc/dXa =2*Dense(Mu,nu)*(Sum over mu centered on A)(Sumover all nu)
  !  Integral((df/drhoa dPhimu/dXA Phinu)-
  !  (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b))
  !  DOT Grad(dPhimu/dXa Phinu))
  !
  !  where F alpha mu nu is the the alpha spin portion of the operator matrix
  !  element mu, nu,
  !  df/drhoa is the derivative of the functional by the alpha density,
  !  df/dgaa is the derivative of the functional by the alpha gradient
  !  invariant, i.e. the dot product of the gradient of the alpha
  !  density with itself.
  !  df/dgab is the derivative of the functional by the dot product of
  !  the gradient of the alpha density with the beta density.
  !  Grad(Phimu Phinu) is the gradient of Phimu times Phinu.
  !-------------------------------------------------------------------------
  
     use allmod
     use quick_gridpoints_module
     use quick_dft_module, only: b3lypf, b3lyp_e, becke, becke_e, lyp, lyp_e
     use xc_f90_types_m
     use xc_f90_lib_m
     implicit none
  
     integer :: iatm, ibas, ibin, icount, ifunc, igp, jbas, jcount, ibasstart, irad_init, &
     irad_end, ierror, imomentum
     double precision :: density, densityb, densitysum, dfdgaa, dfdgaa2, dfdgab, dfdgbb, &
     dfdgab2, dfdr, dfdrb, dfdr2, dphi2dx, dphi2dy, dphi2dz, dphidx, dphidy, dphidz, &
     gax, gay, gaz, gbx, gby, gbz, gaa, gab, gbb, gridx, gridy, gridz, phi, phi2, quicktest, &
     sigma, sswt, temp, tempgx, tempgy, tempgz, tsttmp_exc, tsttmp_vrhoa, &
     tsttmp_vsigmaa, weight, xdot, ydot, zdot, xiaodot, zkec, Ex, Ec, Eelxc, excpp, &
     xdotb, ydotb, zdotb, dxdx, dxdy, dxdz, dydy, dydz, dzdz
     double precision, dimension(2) :: libxc_rho
     double precision, dimension(3) :: libxc_sigma
     double precision, dimension(1) :: libxc_exc
     double precision, dimension(2) :: libxc_vrho
     double precision, dimension(3) :: libxc_vsigma
     type(xc_f90_pointer_t), dimension(quick_method%nof_functionals) ::xc_func
     type(xc_f90_pointer_t), dimension(quick_method%nof_functionals) ::xc_info
  
#ifdef MPIV
     include "mpif.h"
#endif
  
#if defined CUDA || defined CUDA_MPIV
  
     if(quick_method%bCUDA) then
  
        call gpu_reupload_dft_grid()
  
#ifdef OSHELL
        call gpu_get_oshell_xcgrad(quick_qm_struct%gradient)
#else
        call gpu_get_cshell_xcgrad(quick_qm_struct%gradient)
#endif
  
        call gpu_delete_dft_grid()
  
     endif
#else
  
     if(quick_method%uselibxc) then
  !  Initiate the libxc functionals
        do ifunc=1, quick_method%nof_functionals
#ifdef OSHELL
              call xc_f90_func_init(xc_func(ifunc), xc_info(ifunc), &
              quick_method%functional_id(ifunc),XC_POLARIZED)
#else
              call xc_f90_func_init(xc_func(ifunc), &
              xc_info(ifunc),quick_method%functional_id(ifunc),XC_UNPOLARIZED)
#endif
        enddo
     endif
  
#if defined MPIV && !defined CUDA_MPIV
        if(bMPI) then
           irad_init = quick_dft_grid%igridptll(mpirank+1)
           irad_end = quick_dft_grid%igridptul(mpirank+1)
        else
           irad_init = 1
           irad_end = quick_dft_grid%nbins
        endif
        do Ibin=irad_init, irad_end
#else
        do Ibin=1, quick_dft_grid%nbins
#endif
  
  !  Calculate the weight of the grid point in the SSW scheme.  If
  !  the grid point has a zero weight, we can skip it.
  
          Igp=quick_dft_grid%bin_counter(Ibin)+1
  
          do while(Igp < quick_dft_grid%bin_counter(Ibin+1)+1)
  
             gridx=quick_dft_grid%gridxb(Igp)
             gridy=quick_dft_grid%gridyb(Igp)
             gridz=quick_dft_grid%gridzb(Igp)
  
             sswt=quick_dft_grid%gridb_sswt(Igp)
             weight=quick_dft_grid%gridb_weight(Igp)
             Iatm=quick_dft_grid%gridb_atm(Igp)
              
              if (weight < quick_method%DMCutoff ) then
                 continue
              else
  
                 icount=quick_dft_grid%basf_counter(Ibin)+1
                 do while (icount < quick_dft_grid%basf_counter(Ibin+1)+1)
                    Ibas=quick_dft_grid%basf(icount)+1
  
                    call pteval_new_imp(gridx,gridy,gridz,phi,dphidx,dphidy, &
                    dphidz,Ibas,icount)
  
                    phixiao(Ibas)=phi
                    dphidxxiao(Ibas)=dphidx
                    dphidyxiao(Ibas)=dphidy
                    dphidzxiao(Ibas)=dphidz
  
                    icount=icount+1
                 enddo
  
                 
  
  !  evaluate the densities at the grid point and the gradient at that grid point            
#ifdef OSHELL
                 call denspt_oshell(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                 gbx,gby,gbz,Ibin)
#else          
                 call denspt_cshell(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                 gbx,gby,gbz,Ibin)
#endif
  
#ifdef OSHELL
                 if ((density < quick_method%DMCutoff) .and. (densityb < quick_method%DMCutoff)) then
#else
                 if (density < quick_method%DMCutoff ) then
#endif
                    continue
  
                 else
  !  This allows the calculation of the derivative of the functional
  !  with regard to the density (dfdr), with regard to the alpha-alpha
  !  density invariant (df/dgaa), and the alpha-beta density invariant.
  
#ifdef OSHELL
                    gaa = (gax*gax+gay*gay+gaz*gaz)
                    gbb = (gbx*gbx+gby*gby+gbz*gbz)
                    gab = (gax*gbx+gay*gby+gaz*gbz)
        
                    libxc_rho(1)=density
                    libxc_rho(2)=densityb
        
                    libxc_sigma(1)=gaa
                    libxc_sigma(2)=gab
                    libxc_sigma(3)=gbb
#else
                    densitysum=2.0d0*density
                    sigma=4.0d0*(gax*gax+gay*gay+gaz*gaz)
  
                    libxc_rho(1)=densitysum
                    libxc_sigma(1)=sigma
#endif
  
                    excpp=0.0d0
                    dfdr=0.0d0
                    dfdrb=0.0d0
                
                    dfdgaa=0.0d0
                    dfdgab=0.0d0
                    dfdgbb=0.0d0
  
                    if(quick_method%uselibxc) then
                       do ifunc=1, quick_method%nof_functionals
                          select case(xc_f90_info_family(xc_info(ifunc)))
                             case(XC_FAMILY_LDA)
                                call xc_f90_lda_exc_vxc(xc_func(ifunc),1,libxc_rho(1), &
                                libxc_exc(1), libxc_vrho(1))
                                libxc_vsigma(1) = 0.0d0
                                libxc_vsigma(2) = 0.0d0
                                libxc_vsigma(3) = 0.0d0
                             case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                                call xc_f90_gga_exc_vxc(xc_func(ifunc),1,libxc_rho(1), libxc_sigma(1), &
                                libxc_exc(1), libxc_vrho(1), libxc_vsigma(1))
                          end select

                        excpp=excpp+libxc_exc(1)
                        dfdr=dfdr+libxc_vrho(1)
                        dfdgaa=dfdgaa+libxc_vsigma(1)
#ifdef OSHELL
                        dfdrb=dfdrb+libxc_vrho(2)
                        dfdgab=dfdgab+libxc_vsigma(2)
                        dfdgbb=dfdgbb+libxc_vsigma(3)
#endif  

                       enddo
  
                       zkec=(density+densityb)*excpp
#ifdef OSHELL
                       xdot = 2.d0*dfdgaa*gax + dfdgab*gbx
                       ydot = 2.d0*dfdgaa*gay + dfdgab*gby
                       zdot = 2.d0*dfdgaa*gaz + dfdgab*gbz
                    
                       xdotb = 2.d0*dfdgbb*gbx + dfdgab*gax
                       ydotb = 2.d0*dfdgbb*gby + dfdgab*gay
                       zdotb = 2.d0*dfdgbb*gbz + dfdgab*gaz
#else               
                       xdot = 4.0d0*dfdgaa*gax
                       ydot = 4.0d0*dfdgaa*gay
                       zdot = 4.0d0*dfdgaa*gaz
                    
                    elseif(quick_method%BLYP) then
  
                       call becke_E(density, densityb, gax, gay, gaz, gbx, gby,gbz, Ex)
                       call lyp_e(density, densityb, gax, gay, gaz, gbx, gby, gbz,Ec)
  
                       zkec=Ex+Ec
  
                       call becke(density, gax, gay, gaz, gbx, gby, gbz, dfdr, dfdgaa, dfdgab)
                       call lyp(density, densityb, gax, gay, gaz, gbx, gby, gbz, dfdr2, dfdgaa2, dfdgab2)
              
                       dfdr = dfdr + dfdr2
                       dfdgaa = dfdgaa + dfdgaa2
                       dfdgab = dfdgab + dfdgab2
  
                       xdot = 2.d0*dfdgaa*gax + dfdgab*gbx
                       ydot = 2.d0*dfdgaa*gay + dfdgab*gby
                       zdot = 2.d0*dfdgaa*gaz + dfdgab*gbz
  
                    elseif(quick_method%B3LYP) then
  
                       call b3lyp_e(densitysum, sigma, zkec)
                       call b3lypf(densitysum, sigma, dfdr, xiaodot)
  
                       xdot=xiaodot*gax
                       ydot=xiaodot*gay
                       zdot=xiaodot*gaz
#endif 
                    endif
  
  ! Now loop over basis functions and compute the addition to the matrix
  ! element.
                    icount=quick_dft_grid%basf_counter(Ibin)+1
                    do while (icount < quick_dft_grid%basf_counter(Ibin+1)+1)
                       Ibas=quick_dft_grid%basf(icount)+1
  
                       phi=phixiao(Ibas)
                       dphidx=dphidxxiao(Ibas)
                       dphidy=dphidyxiao(Ibas)
                       dphidz=dphidzxiao(Ibas)
  
                       !call pteval_new_imp(gridx,gridy,gridz,phi,dphidx,dphidy, &
                       !dphidz,Ibas,icount)
  
  
                       quicktest = DABS(dphidx+dphidy+dphidz+phi)
                       
                       if (quicktest < quick_method%DMCutoff ) then
                          continue
                       else
                          call pt2der(gridx,gridy,gridz,dxdx,dxdy,dxdz, &
                          dydy,dydz,dzdz,Ibas,icount)
  
                          Ibasstart=(quick_basis%ncenter(Ibas)-1)*3
  
                          jcount=quick_dft_grid%basf_counter(Ibin)+1
                          do while(jcount<quick_dft_grid%basf_counter(Ibin+1)+1)
                             Jbas = quick_dft_grid%basf(jcount)+1 
  
                             phi2=phixiao(Jbas)
                             dphi2dx=dphidxxiao(Jbas)
                             dphi2dy=dphidyxiao(Jbas)
                             dphi2dz=dphidzxiao(Jbas)
  
                             !call pteval_new_imp(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
                             !dphi2dz,Jbas,jcount)
  
                             quick_qm_struct%gradient(Ibasstart+1) =quick_qm_struct%gradient(Ibasstart+1) - &
                             2.d0*quick_qm_struct%dense(Ibas,Jbas)*weight*&
                             (dfdr*dphidx*phi2 &
                             + xdot*(dxdx*phi2+dphidx*dphi2dx) &
                             + ydot*(dxdy*phi2+dphidx*dphi2dy) &
                             + zdot*(dxdz*phi2+dphidx*dphi2dz))
                             quick_qm_struct%gradient(Ibasstart+2)= quick_qm_struct%gradient(Ibasstart+2) - &
                             2.d0*quick_qm_struct%dense(Ibas,Jbas)*weight*&
                             (dfdr*dphidy*phi2 &
                             + xdot*(dxdy*phi2+dphidy*dphi2dx) &
                             + ydot*(dydy*phi2+dphidy*dphi2dy) &
                             + zdot*(dydz*phi2+dphidy*dphi2dz))
                             quick_qm_struct%gradient(Ibasstart+3)= quick_qm_struct%gradient(Ibasstart+3) - &
                             2.d0*quick_qm_struct%dense(Ibas,Jbas)*weight*&
                             (dfdr*dphidz*phi2 &
                             + xdot*(dxdz*phi2+dphidz*dphi2dx) &
                             + ydot*(dydz*phi2+dphidz*dphi2dy) &
                             + zdot*(dzdz*phi2+dphidz*dphi2dz))

#ifdef OSHELL
                             quick_qm_struct%gradient(Ibasstart+1) =quick_qm_struct%gradient(Ibasstart+1) - &
                             2.d0*quick_qm_struct%denseb(Ibas,Jbas)*weight*&
                             (dfdrb*dphidx*phi2 &
                             + xdotb*(dxdx*phi2+dphidx*dphi2dx) &
                             + ydotb*(dxdy*phi2+dphidx*dphi2dy) &
                             + zdotb*(dxdz*phi2+dphidx*dphi2dz))
                             quick_qm_struct%gradient(Ibasstart+2)= quick_qm_struct%gradient(Ibasstart+2) - &
                             2.d0*quick_qm_struct%denseb(Ibas,Jbas)*weight*&
                             (dfdrb*dphidy*phi2 &
                             + xdotb*(dxdy*phi2+dphidy*dphi2dx) &
                             + ydotb*(dydy*phi2+dphidy*dphi2dy) &
                             + zdotb*(dydz*phi2+dphidy*dphi2dz))
                             quick_qm_struct%gradient(Ibasstart+3)= quick_qm_struct%gradient(Ibasstart+3) - &
                             2.d0*quick_qm_struct%denseb(Ibas,Jbas)*weight*&
                             (dfdrb*dphidz*phi2 &
                             + xdotb*(dxdz*phi2+dphidz*dphi2dx) &
                             + ydotb*(dydz*phi2+dphidz*dphi2dy) &
                             + zdotb*(dzdz*phi2+dphidz*dphi2dz))
#endif

                             jcount=jcount+1
                          enddo
                       endif
  
                    icount=icount+1
                    enddo
  
  !  We are now completely done with the derivative of the exchange correlation energy with nuclear displacement
  !  at this point. Now we need to do the quadrature weight derivatives. At this point in the loop, we know that
  !  the density and the weight are not zero. Now check to see fi the weight is one. If it isn't, we need to
  !  actually calculate the energy and the derivatives of the quadrature at this point. Due to the volume of code,
  !  this is done in sswder. Note that if a new weighting scheme is ever added, this needs
  !  to be modified with a second subprogram.
                    if (sswt == 1.d0) then
                       continue
                    else
                       call sswder(gridx,gridy,gridz,zkec,weight/sswt,Iatm)
                    endif
                 endif
              endif
  !         enddo
  
        Igp=Igp+1
        enddo
     enddo
  
     if(quick_method%uselibxc) then
  !  Uninitilize libxc functionals
        do ifunc=1, quick_method%nof_functionals
           call xc_f90_func_end(xc_func(ifunc))
        enddo
     endif
#endif

!#ifdef DEBUG
!    if (quick_method%debug) then
          write (*,'(/," DEBUG STEP 5 :  XC GRADIENT ADDED:")')
          do Iatm=1,natom
              do Imomentum=1,3
                  write (*,'(I5,7x,F20.10)')Iatm, &
                  quick_qm_struct%gradient((Iatm-1)*3+Imomentum)
              enddo
          enddo
!    endif
!#endif  

     return
  
#ifdef OSHELL
  end subroutine get_oshell_xc_grad
#else
  end subroutine get_cshell_xc_grad
#endif
  

#ifndef OSHELL  
  subroutine get_ijbas_derivative(Imomentum, Ibas, Jbas, mbas, mstart, ijcon, DENSEJI)
  
  !-------------------------------------------------------------------------
  !  The purpose of this subroutine is to compute the I and J basis function
  !  derivatives required for get_kinetic_grad subroutine. The input variables
  !  mbas, mstart and ijcon are used to differentiate between I and J
  !  basis functions. For I basis functions, ijcon should be  true and should
  !  be false for J.  
  !-------------------------------------------------------------------------   
     use allmod
     use quick_overlap_module, only: gpt, opf, overlap
     use quick_oei_module, only: ekinetic
     implicit double precision(a-h,o-z)
     logical :: ijcon   
     double precision g_table(200), valopf
     integer i,j,k,ii,jj,kk,g_count
  
     dSM = 0.0d0
     dKEM = 0.0d0
  
     Ax = xyz(1,quick_basis%ncenter(Jbas))
     Bx = xyz(1,quick_basis%ncenter(Ibas))
     Ay = xyz(2,quick_basis%ncenter(Jbas))
     By = xyz(2,quick_basis%ncenter(Ibas))
     Az = xyz(3,quick_basis%ncenter(Jbas))
     Bz = xyz(3,quick_basis%ncenter(Ibas))
     
     itype(Imomentum,mbas) = itype(Imomentum,mbas)+1
  
     ii = itype(1,Ibas)
     jj = itype(2,Ibas)
     kk = itype(3,Ibas)
     i = itype(1,Jbas)
     j = itype(2,Jbas)
     k = itype(3,Jbas)
     g_count = i+ii+j+jj+k+kk+2
  
     do Icon=1,ncontract(Ibas)
        b = aexp(Icon,Ibas)
        do Jcon=1,ncontract(Jbas)
           a = aexp(Jcon,Jbas)
  
           valopf = opf(a, b, dcoeff(Jcon,Jbas), dcoeff(Icon,Ibas), Ax,&
                    Ay, Az, Bx, By, Bz)
  
           if(abs(valopf) .gt. quick_method%coreIntegralCutoff) then                  
  
             call gpt(a,b,Ax,Ay,Az,Bx,By,Bz,Px,Py,Pz,g_count,g_table)      
            
             if(ijcon) then
                mcon = Icon
             else
                mcon = Jcon
             endif
             dSM= dSM + 2.d0*aexp(mcon,mbas)* &
             dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
             *overlap(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az,Bx,By,Bz,Px,Py,Pz,g_table)
             dKEM = dKEM + 2.d0*aexp(mcon,mbas)* &
             dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
             *ekinetic(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az,Bx,By,Bz,Px,Py,Pz,g_table)
           endif
        enddo
     enddo
  
     itype(Imomentum,mbas) = itype(Imomentum,mbas)-1
  
     if (itype(Imomentum,mbas) /= 0) then
        itype(Imomentum,mbas) = itype(Imomentum,mbas)-1
  
        ii = itype(1,Ibas)
        jj = itype(2,Ibas)
        kk = itype(3,Ibas)
        i = itype(1,Jbas)
        j = itype(2,Jbas)
        k = itype(3,Jbas)
        g_count = i+ii+j+jj+k+kk+2
  
        do Icon=1,ncontract(Ibas)
           b = aexp(Icon,Ibas)
           do Jcon=1,ncontract(Jbas)
             a = aexp(Jcon,Jbas)
             
             valopf = opf(a, b, dcoeff(Jcon,Jbas), dcoeff(Icon,Ibas), Ax,&
                      Ay, Az, Bx, By, Bz)
             
             if(abs(valopf) .gt. quick_method%coreIntegralCutoff) then
             
               call gpt(a,b,Ax,Ay,Az,Bx,By,Bz,Px,Py,Pz,g_count,g_table)      
             
               dSM = dSM - dble(itype(Imomentum,mbas)+1)* &
               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
               *overlap(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az,Bx,By,Bz,Px,Py,Pz,g_table)
               dKEM = dKEM - dble(itype(Imomentum,mbas)+1)* &
               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
               *ekinetic(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az,Bx,By,Bz,Px,Py,Pz,g_table)
             endif
           enddo
        enddo
  
        itype(Imomentum,mbas) = itype(Imomentum,mbas)+1
     endif
  
     quick_qm_struct%gradient(mstart+Imomentum) = quick_qm_struct%gradient(mstart+Imomentum) &
     -dSM*quick_scratch%hold(Jbas,Ibas)*2.d0 &
     +dKEM*DENSEJI*2.d0
  
     return
  
  end subroutine get_ijbas_derivative
#endif

#ifdef OSHELL
end module quick_oshell_gradient_module
#else
end module quick_cshell_gradient_module
#endif
