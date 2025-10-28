! subroutine xcrank_coupled

!     use V_proce
!     use V_cons
!     use V_para
!     use V_hom
!     use V_spe
!     use V_crank

!     implicit none

!     integer :: SeBra1, SeBra2, SeKet1, SeKet2
!     !serial number of coupled basis in one certain cranking basis
!     real :: Nb, Lb, Jb, Mb
!     real :: Nk, Lk, Jk, Mk

!     integer :: DIM
!     integer :: EV=1, IROT
!     integer, dimension(1000,1000) :: B,Z

!     real, allocatable, dimension(:,:) :: Hr !H_real
!     real, allocatable, dimension(:,:) :: Hi !H_imaginary

!     real, allocatable, dimension(:) :: Etmp
!     real, allocatable, dimension(:,:) :: Htmp
!     real, allocatable, dimension(:,:) :: Stmp

!     DIM=CrankDim(N)

!     allocate(Hr(DIM,DIM))

!     signature=-0.5

!     do bra=1, CrankDim(N)
!     do ket=bra, CrankDim(N)

!         SeBra1=int(CrankBV(N,bra,1))
!         SeBra2=int(CrankBV(N,bra,2))
!         SeKet1=int(CrankBV(N,ket,1))
!         SeKet2=int(CrankBV(N,ket,2))

!         Nb=CBasis(N,SeBra1,1)
!         Lb=CBasis(N,SeBra1,2)
!         Jb=CBasis(N,SeBra1,3)
!         Mb=CBasis(N,SeBra1,4)

!         Nk=CBasis(N,SeKet1,1)
!         Lk=CBasis(N,SeKet1,2)
!         Jk=CBasis(N,SeKet1,3)
!         Mk=CBasis(N,SeKet1,4)

!         ! write(*,'(2i4,100f6.1)') bra,ket, Nb,Nk, Lb,Lk, Jb,Jk, Mb,Mk

!         Hr(bra,ket)=Hspc(N,SeBra1,SeKet1) &
!         + ((-1)**(Jb+signature))*Hspc(N,SeBra2,SeKet1) &
!         + ((-1)**(Jk+signature))*Hspc(N,SeBra1,SeKet2) &
!         + ((-1)**(Jb+Jk+2*signature))*Hspc(N,SeBra2,SeKet2)

!         Hr(bra,ket)=0.5*Hr(bra,ket)

!         ! write(*,'(100f12.6)') Hspc(N,SeBra1,SeKet1), (-1)**(Jb+Jk), Hspc(N,SeBra2,SeKet2)

!         if (Nb/=Nk .or. Lb/=Lk .or. Jb/=Jk ) cycle

!         if (Mb==Mk+1.) then
!             Hr(bra,ket)=Hr(bra,ket)-omega*0.25*sqrt(Jk*(Jk+1)-Mk*(Mk+1))
!             ! if (bra==1 .and. ket==2) 
!             ! write(*,'(3I4,2f6.1, 2f12.6)') 1, bra,ket, Jk,Mk, sqrt(Jk*(Jk+1)-Mk*(Mk+1)), Hcrk(N,bra,ket)
!         end if

!         if (Mb==Mk-1.) then
!             Hr(bra,ket)=Hr(bra,ket)-omega*0.25*sqrt(Jk*(Jk+1)-Mk*(Mk-1))
!             ! if (bra==1 .and. ket==2) 
!             ! write(*,'(3I4,2f6.1, 2f12.6)') 2, bra,ket, Jk,Mk, sqrt(Jk*(Jk+1)-Mk*(Mk-1)), Hcrk(N,bra,ket)
!         end if

!         if (-Mb==Mk+1.) then
!             Hr(bra,ket)=Hr(bra,ket)-omega*0.25*((-1)**(Jb+signature)) &
!             *sqrt(Jk*(Jk+1)-Mk*(Mk+1))
!             ! if (bra==1 .and. ket==2) 
!             ! write(*,'(3I4,2f6.1, 2f12.6)') 3, bra,ket, Jk,Mk, sqrt(Jk*(Jk+1)-Mk*(Mk+1)), Hcrk(N,bra,ket)
!         end if

!         if (-Mb==Mk-1.) then
!             Hr(bra,ket)=Hr(bra,ket)-omega*0.25*((-1)**(Jb+signature)) &
!             *sqrt(Jk*(Jk+1)-Mk*(Mk-1))
!             ! if (bra==1 .and. ket==2) 
!             ! write(*,'(3I4,2f6.1, 2f12.6)') 4, bra,ket, Jk,Mk, sqrt(Jk*(Jk+1)-Mk*(Mk-1)), Hcrk(N,bra,ket)
!         end if

!         if (Mb==-Mk+1.) then
!             Hr(bra,ket)=Hr(bra,ket)-omega*0.25*((-1)**(Jk+signature)) &
!             *sqrt(Jk*(Jk+1)+Mk*(-Mk+1))
!             ! if (bra==1 .and. ket==2) 
!             ! write(*,'(3I4,2f6.1, 2f12.6)') 5, bra,ket, Jk,Mk, sqrt(Jk*(Jk+1)+Mk*(-Mk+1)), Hcrk(N,bra,ket)
!         end if

!         if (Mb==-Mk-1.) then
!             Hr(bra,ket)=Hr(bra,ket)-omega*0.25*((-1)**(Jk+signature)) &
!             *sqrt(Jk*(Jk+1)+Mk*(-Mk-1))
!             ! if (bra==1 .and. ket==2) 
!             ! write(*,'(3I4,2f6.1, 2f12.6)') 6, bra,ket, Jk,Mk, sqrt(Jk*(Jk+1)+Mk*(-Mk-1)), Hcrk(N,bra,ket)
!         end if

!         if (-Mb==-Mk+1.) then
!             Hr(bra,ket)=Hr(bra,ket)-omega*0.25*((-1)**(Jb+Jk+2*signature)) &
!             *sqrt(Jk*(Jk+1)+Mk*(-Mk+1))
!             ! if (bra==1 .and. ket==2) 
!             ! write(*,'(3I4,2f6.1, 2f12.6)') 7, bra,ket, Jk,Mk, sqrt(Jk*(Jk+1)+Mk*(-Mk+1)), Hcrk(N,bra,ket)
!         end if

!         if (-Mb==-Mk-1.) then
!             Hr(bra,ket)=Hr(bra,ket)-omega*0.25*((-1)**(Jb+Jk+2*signature)) &
!             *sqrt(Jk*(Jk+1)+Mk*(-Mk-1))
!             ! if (bra==1 .and. ket==2) 
!             ! write(*,'(3I4,2f6.1, 2f12.6)') 8, bra,ket, Jk,Mk, sqrt(Jk*(Jk+1)+Mk*(-Mk-1)), Hcrk(N,bra,ket)
!         end if

!         Hr(ket,bra)=Hr(bra,ket)
!     end do
!     end do

!     DIM=CrankDim(N) 

!     ! do bra=1, DIM
!     !     write(*,'(1000f12.6)') Hr(bra,:)
!     ! end do

!     allocate(Htmp(DIM,DIM))
!     allocate(Etmp(DIM))
!     allocate(Stmp(DIM,DIM))

!     call diagonalise(DIM,DIM,EV,Hr,Etmp,Stmp,IROT,B,Z)

!     ! write(*,*) 'signature=', signature
!     ! do count=1, DIM
!     !     write(*,*) Etmp(count)/hbom0
!     ! end do

!     do count=1, CrankDim(N)
!         index=index+1
!         CrkE(index)=real(Etmp(count))
!         CrkBasis(index,:)=CBasis(N,int(CrankBV(N,count,1)),:)
!         CrkWF(:,index)=Stmp(:,count)
!     end do

!     deallocate(Stmp)
!     deallocate(Etmp)
!     deallocate(Htmp)

!     deallocate(Hr)

! end subroutine xcrank_coupled