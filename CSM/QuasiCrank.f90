! subroutine QuasiCrank

!     use V_proce
!     use V_cons
!     use V_para
!     use V_hom
!     use V_spe
!     use V_quasi
!     use V_crank

!     implicit none

!     integer :: DIM
!     integer :: braI, ketI
!     integer :: SeBra1, SeBra2, SeKet1, SeKet2
!     real :: Nb, Lb, Jb, Mb
!     real :: Nk, Lk, Jk, Mk

!     integer :: EV=1, IROT
!     integer, dimension(1000,1000) :: B,Z

!     real, allocatable, dimension(:) :: Etmp !single particle energy
!     real, allocatable, dimension(:,:) :: Htmp !single particle hamiltonian
!     real, allocatable, dimension(:,:) :: Stmp !single particle wave function

!     DIM=CBasisDim(N)/2

!     Efermi=6.485*hbom0
    
!     if (omega/hbom0 < 0.05) then
!         Delta=0.12*hbom0
!     else
!         Delta=-2.4*(omega/hbom0)+0.24
!         Delta=Delta*hbom0
!     end if

!     ! Efermi=0
!     ! Delta=0

!     do bra=1, DIM
!         braI=bra+DIM
!     do ket=1, DIM
!         ketI=ket+DIM

!         SeBra1=int(CrankBV(N,bra,1))
!         SeBra2=int(CrankBV(N,bra,2))
!         SeKet1=int(CrankBV(N,ket,1))
!         SeKet2=int(CrankBV(N,ket,2))

!         write(*,'(10f4.1)') CBasis(N,SeBra1,4), CBasis(N,SeBra2,4), &
!         CBasis(N,SeKet1,4), CBasis(N,SeKet2,4)

!         Nb=CBasis(N,SeBra1,1)
!         Lb=CBasis(N,SeBra1,2)
!         Jb=CBasis(N,SeBra1,3)
!         Mb=CBasis(N,SeBra1,4)

!         Nk=CBasis(N,SeKet1,1)
!         Lk=CBasis(N,SeKet1,2)
!         Jk=CBasis(N,SeKet1,3)
!         Mk=CBasis(N,SeKet1,4)
        
!         Hcrk(N,bra,ket)=Hspc(N,SeBra1,SeKet1) &
!         +((-1)**(Jk+0.5))*Hspc(N,SeBra1,SeKet2) &
!         +((-1)**(Jb+0.5))*Hspc(N,SeBra2,SeKet1) &
!         +((-1)**(Jb+Jk+1))*Hspc(N,SeBra2,SeKet2)
!         Hcrk(N,bra,ket)=0.5*Hcrk(N,bra,ket)

!         ! Hcrk(N,braI,ketI)=Hspc(N,SeBra1,SeKet1) &
!         ! +((-1)**(Jk-0.5))*Hspc(N,SeBra1,SeKet2) &
!         ! +((-1)**(Jb-0.5))*Hspc(N,SeBra2,SeKet1) &
!         ! +((-1)**(Jb+Jk-1))*Hspc(N,SeBra2,SeKet2)
!         ! Hcrk(N,braI,ketI)=0.5*Hcrk(N,braI,ketI)
!         Hcrk(N,braI,ketI)=-Hcrk(N,bra,ket)

!         if (Nb/=Nk .or. Lb/=Lk .or. Jb/=Jk) cycle 

!         !!!\alpha=1/2!!!
!         if (Mb==Mk+1) &
!         Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
!         -0.25*sqrt(Jk*(Jk+1)-Mk*(Mk+1))*omega
!         if (Mb==Mk-1) &
!         Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
!         -0.25*sqrt(Jk*(Jk+1)-Mk*(Mk-1))*omega

!         if (Mb==-Mk+1) &
!         Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
!         -0.25*((-1)**(Jk+0.5))*sqrt(Jk*(Jk+1)+Mk*(-Mk+1))*omega
!         if (Mb==-Mk-1) &
!         Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
!         -0.25*((-1)**(Jk+0.5))*sqrt(Jk*(Jk+1)+Mk*(-Mk-1))*omega

!         if (-Mb==Mk+1) &
!         Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
!         -0.25*((-1)**(Jb+0.5))*sqrt(Jk*(Jk+1)-Mk*(Mk+1))*omega
!         if (-Mb==Mk-1) &
!         Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
!         -0.25*((-1)**(Jb+0.5))*sqrt(Jk*(Jk+1)-Mk*(Mk-1))*omega

!         if (-Mb==-Mk+1) &
!         Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
!         -0.25*((-1)**(Jb+Jk+1))*sqrt(Jk*(Jk+1)+Mk*(-Mk+1))*omega
!         if (-Mb==-Mk-1) &
!         Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
!         -0.25*((-1)**(Jb+Jk+1))*sqrt(Jk*(Jk+1)+Mk*(-Mk-1))*omega


!         !!!\alpha=-1/2!!!
!         ! if (Mb==Mk+1) &
!         ! Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         ! -0.25*sqrt(Jk*(Jk+1)-Mk*(Mk+1))*omega
!         ! if (Mb==Mk-1) &
!         ! Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         ! -0.25*sqrt(Jk*(Jk+1)-Mk*(Mk-1))*omega

!         ! if (Mb==-Mk+1) &
!         ! Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         ! -0.25*((-1)**(Jk-0.5))*sqrt(Jk*(Jk+1)+Mk*(-Mk+1))*omega
!         ! if (Mb==-Mk-1) &
!         ! Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         ! -0.25*((-1)**(Jk-0.5))*sqrt(Jk*(Jk+1)+Mk*(-Mk-1))*omega

!         ! if (-Mb==Mk+1) &
!         ! Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         ! -0.25*((-1)**(Jb-0.5))*sqrt(Jk*(Jk+1)-Mk*(Mk+1))*omega
!         ! if (-Mb==Mk-1) &
!         ! Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         ! -0.25*((-1)**(Jb-0.5))*sqrt(Jk*(Jk+1)-Mk*(Mk-1))*omega

!         ! if (-Mb==-Mk+1) &
!         ! Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         ! -0.25*((-1)**(Jb+Jk-1))*sqrt(Jk*(Jk+1)+Mk*(-Mk+1))*omega
!         ! if (-Mb==-Mk-1) &
!         ! Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         ! -0.25*((-1)**(Jb+Jk-1))*sqrt(Jk*(Jk+1)+Mk*(-Mk-1))*omega

!         if (Mb==Mk+1) &
!         Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         +0.25*sqrt(Jk*(Jk+1)-Mk*(Mk+1))*omega
!         if (Mb==Mk-1) &
!         Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         +0.25*sqrt(Jk*(Jk+1)-Mk*(Mk-1))*omega

!         if (Mb==-Mk+1) &
!         Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         +0.25*((-1)**(Jk+0.5))*sqrt(Jk*(Jk+1)+Mk*(-Mk+1))*omega
!         if (Mb==-Mk-1) &
!         Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         +0.25*((-1)**(Jk+0.5))*sqrt(Jk*(Jk+1)+Mk*(-Mk-1))*omega

!         if (-Mb==Mk+1) &
!         Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         +0.25*((-1)**(Jb+0.5))*sqrt(Jk*(Jk+1)-Mk*(Mk+1))*omega
!         if (-Mb==Mk-1) &
!         Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         +0.25*((-1)**(Jb+0.5))*sqrt(Jk*(Jk+1)-Mk*(Mk-1))*omega

!         if (-Mb==-Mk+1) &
!         Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         +0.25*((-1)**(Jb+Jk+1))*sqrt(Jk*(Jk+1)+Mk*(-Mk+1))*omega
!         if (-Mb==-Mk-1) &
!         Hcrk(N,braI,ketI)=Hcrk(N,braI,ketI) &
!         +0.25*((-1)**(Jb+Jk+1))*sqrt(Jk*(Jk+1)+Mk*(-Mk-1))*omega

!     end do  
!         Hcrk(N,bra,bra)=Hcrk(N,bra,bra)-Efermi
!         Hcrk(N,braI,braI)=Hcrk(N,braI,braI)+Efermi
!         Hcrk(N,bra,braI)=Delta
!         Hcrk(N,braI,bra)=-Delta
!     end do

!     stop

!     DIM=DIM*2
!     allocate(Htmp(DIM,DIM))
!     allocate(Etmp(DIM))
!     allocate(Stmp(DIM,DIM))

!     Htmp=0
!     do bra=1, DIM
!     do ket=1, DIM
!         Htmp(bra,ket)=Hcrk(N,bra,ket)
!     end do
!     end do

!     ! do bra=1, DIM
!     !     write(*,'(1000f12.6)') Htmp(bra,:)
!     ! end do

!     ! stop

!     call diagonalise(DIM,DIM,EV,Htmp,Etmp,Stmp,IROT,B,Z)

!     do count=1, DIM
!         index=index+1
!         CrkE(index)=Etmp(count)
!         CrkWF(:,index)=Stmp(:,index)
!     end do

!     ! do count=1, DIM
!     !     write(*,*) Etmp(count)/hbom0
!     ! end do

!     deallocate(Stmp)
!     deallocate(Etmp)
!     deallocate(Htmp)

! end subroutine QuasiCrank