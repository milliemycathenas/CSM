subroutine NilSPE
    !Nil_spe calculates the single particle energies in the Nilsson potential
    
    use V_proce
    use V_cons
    use V_para
    use V_hom
    use V_spe
    use V_crank

    implicit none

    integer :: DIM
    integer :: EV=1, IROT
    integer, dimension(1000,1000) :: B,Z

    !!!variables for CG factor subroutine start!!!
        double precision :: J1, J2, J3, M1, M2, M3
        double precision :: CG200,CG2mm,CG400,CG4mm
        double precision :: J1bra, J1ket, J3bra, J3ket
        double precision :: M5bra, M5ket, M6bra, M6ket, M3bra, M3ket 
        !J5 for \Sigma=1/2, J6 for \Sigma=-1/2
        double precision :: CGbra5, CGbra6, CGket5, CGket6
    !!!variables for CG factor subroutine end!!!

    !!!for sSYEV start!!!
        integer(kind=8), parameter :: LWMAX=10000
        integer(kind=8) :: INFO, LWORK
        real(kind=4), dimension(LWMAX) :: work
    !!!for sSYEV end!!!

    !!!process variables start!!!
        integer :: SerBra5, SerBra6, SerKet5, SerKet6
        !index number of the combination of uncoupled basis to coupled basis

        real :: Nb,Lb,Lamb,Sigb, Nk,Lk,Lamk,Sigk
        real :: sqll !(-2/3)*\sqrt{\frac{2l+1}{2l'+1}}

        real, allocatable, dimension(:) :: Etmp !single particle energy
        real, allocatable, dimension(:,:) :: Htmp !single particle hamiltonian
        real, allocatable, dimension(:,:) :: Stmp !single particle wave function
    !!!process variables end!!!

    if (PorN==1) then
        hbom0t=hbom0tp
        hbom0=hbom0p
        kappa=kappap(N+1)
        mu=mup(N+1)
    else if (PorN==0) then
        hbom0t=hbom0tn
        hbom0=hbom0n
        kappa=kappan(N+1)
        mu=mun(N+1)
    end if

    ! kappa=0.075
    ! mu=0.4

    do bra=1, UBasisDim(N)
    do ket=bra, UBasisDim(N)

        Nb=UBasis(N,bra,1)
        Lb=UBasis(N,bra,2)
        Lamb=UBasis(N,bra,3)
        Sigb=UBasis(N,bra,4)

        Nk=UBasis(N,ket,1)
        Lk=UBasis(N,ket,2)
        Lamk=UBasis(N,ket,3)
        Sigk=UBasis(N,ket,4)

        if (Lb/=Lk .and. Lb/=Lk+2 .and. Lb/=Lk-2) cycle

        J2=dble(UBasis(N,ket,2))
        J3=dble(UBasis(N,bra,2))
        M2=dble(UBasis(N,ket,3))
        M3=dble(UBasis(N,bra,3))

        sqll=sqrt((2*Lk+1)/(2*Lb+1))

        if (Sigb==Sigk) then
            call CGFactor(dble(2),J2,J3,dble(0),dble(0),dble(0), CG200)
            call CGFactor(dble(4),J2,J3,dble(0),dble(0),dble(0), CG400)
            if (Lamb==Lamk) then
                call CGFactor(dble(2),J2,J3,dble(0),M2,M3, CG2mm)
                call CGFactor(dble(4),J2,J3,dble(0),M2,M3, CG4mm)
                Hspu(N,bra,ket)=hbom0*sqll &
                *((-2./3.)*eps2*cos(gamrad)*CG200*CG2mm+eps4*d40*CG400*CG4mm)
            end if
            if (Lamb==Lamk+2) then
                call CGFactor(dble(2),J2,J3,dble(2),M2,M3, CG2mm)
                call CGFactor(dble(4),J2,J3,dble(2),M2,M3, CG4mm)
                Hspu(N,bra,ket)=hbom0*sqll &
                *((2./3.)*eps2*sin(gamrad)*CG200*CG2mm/sqrt(2.)+eps4*d42*CG400*CG4mm)
            end if
            if (Lamb==Lamk-2) then
                call CGFactor(dble(2),J2,J3,dble(-2),M2,M3, CG2mm)
                call CGFactor(dble(4),J2,J3,dble(-2),M2,M3, CG4mm)
                Hspu(N,bra,ket)=hbom0*sqll &
                *((2./3.)*eps2*sin(gamrad)*CG200*CG2mm/sqrt(2.)+eps4*d42*CG400*CG4mm)
            end if
            if (Lamb==Lamk+4) then
                call CGFactor(dble(4),J2,J3,dble(4),M2,M3, CG4mm)
                Hspu(N,bra,ket)=hbom0*sqll*eps4*d44*CG400*CG4mm
            end if
            if (Lamb==Lamk-4) then
                call CGFactor(dble(4),J2,J3,dble(-4),M2,M3, CG4mm)
                Hspu(N,bra,ket)=hbom0*sqll*eps2*d44*CG400*CG4mm
            end if

            if (Lb==Lk) then
                Hspu(N,bra,ket)=Hspu(N,bra,ket)*(Nk+1.5)
            else if (Lb==Lk-2) then
                Hspu(N,bra,ket)=Hspu(N,bra,ket)*sqrt(Nk-Lk+2.)*sqrt(Nk+Lk+1.)
            else if (Lb==Lk+2) then
                Hspu(N,bra,ket)=Hspu(N,bra,ket)*sqrt(Nk-Lk)*sqrt(Nk+Lk+3.)
            end if
        end if

        if (Sigb==Sigk .and. Lamb==Lamk .and. Lb==Lk) then
            Hspu(N,bra,ket)=Hspu(N,bra,ket)+hbom0*(Nk+1.5) &
            - 2.*kappa*hbom0t*Lamk*Sigk &
            - kappa*mu*hbom0t &
            *(Lk*(Lk+1.)-0.5*Nk*(Nk+3.))
        end if

        if (Sigb==Sigk+1 .and. Lamb==Lamk-1 .and. Lb==Lk) then
            Hspu(N,bra,ket)=-kappa*hbom0t &
            *sqrt(Lk*(Lk+1)-Lamk*(Lamk-1)) &
            *sqrt(0.5*1.5-Sigk*(Sigk+1))
        end if

        if (Sigb==Sigk-1 .and. Lamb==Lamk+1 .and. Lb==Lk) then
            Hspu(N,bra,ket)=-kappa*hbom0t &
            *sqrt(Lk*(Lk+1)-Lamk*(Lamk+1)) &
            *sqrt(0.5*1.5-Sigk*(Sigk-1))
        end if

        Hspu(N,ket,bra)=Hspu(N,bra,ket)

    end do
    end do

    ! INFO=-1
    ! LWORK=-1
    ! call sSYEV('Vectors', 'Upper', DIM, Hspu, DIM, Etmp, WORK, LWORK, INFO)
    ! LWORK = min(LWMAX, int(work(1)))
    ! !solve eigen problem!
    ! call sSYEV('Vectors', 'Upper', DIM, Hspu, DIM, Etmp, WORK, LWORK, INFO)
    ! !check for convergence!
    ! if (INFO.GT.0) then
    ! write(*,*) 'the alogrithm failed to compute eigenvalues.'
    ! stop
    ! end if

    ! DIM=UBasisDim(N)
    ! allocate(Htmp(DIM,DIM))
    ! allocate(Etmp(DIM))
    ! allocate(Stmp(DIM,DIM))

    ! do bra=1, DIM
    ! do ket=1, DIM
    !     Htmp(bra,ket)=Hspu(N,bra,ket)
    ! end do
    ! end do

    ! call diagonalise(DIM,DIM,EV,Htmp,Etmp,Stmp,IROT,B,Z)

    ! do count=1, DIM
    !     index=index+1
    !     SPE(index)=Etmp(count)
    !     SPBasis(index,:)=UBasis(N,count,:)
    !     SPWF(:,index)=SPWF(:,count)
    ! end do

    ! deallocate(Stmp)
    ! deallocate(Etmp)
    ! deallocate(Htmp)

!!! turn uncoupled Nilsson basis into coupled Nilsson basis start !!!
    do bra=1, CBasisDim(N)
    do ket=1, CBasisDim(N)

        J1bra=dble(CBasis(N,bra,2))
        J1ket=dble(CBasis(N,ket,2))
        M5bra=dble(UBasis(N,int(CBasis(N,bra,5)),3))
        M5ket=dble(UBasis(N,int(CBasis(N,ket,5)),3))
        M6bra=dble(UBasis(N,int(CBasis(N,bra,6)),3))
        M6ket=dble(UBasis(N,int(CBasis(N,ket,6)),3))
        J3bra=dble(CBasis(N,bra,3))
        J3ket=dble(CBasis(N,ket,3))
        M3bra=dble(CBasis(N,bra,4))
        M3ket=dble(CBasis(N,ket,4))

        SerBra5=int(CBasis(N,bra,5))
        SerBra6=int(CBasis(N,bra,6))
        SerKet5=int(CBasis(N,ket,5))
        SerKet6=int(CBasis(N,ket,6))

        call CGFactor(J1bra,dble(0.5),J3bra,M5bra,dble(0.5),M3bra, CGbra5)
        call CGFactor(J1bra,dble(0.5),J3bra,M6bra,dble(-0.5),M3bra, CGbra6)
        call CGFactor(J1ket,dble(0.5),J3ket,M5ket,dble(0.5),M3ket, CGket5)
        call CGFactor(J1ket,dble(0.5),J3ket,M6ket,dble(-0.5),M3ket, CGket6)

        if (SerBra5>0 .and. SerKet5>0) &
        Hspc(N,bra,ket)=CGbra5*CGket5*Hspu(N,SerBra5,SerKet5)
        if (SerBra5>0 .and. SerKet6>0) &
        Hspc(N,bra,ket)=Hspc(N,bra,ket) &
        + CGbra5*CGket6*Hspu(N,SerBra5,SerKet6)
        if (SerBra6>0 .and. SerKet5>0) &
        Hspc(N,bra,ket)=Hspc(N,bra,ket) &
        + CGbra6*CGket5*Hspu(N,SerBra6,SerKet5)
        if (SerBra6>0 .and. SerKet6>0) &
        Hspc(N,bra,ket)=Hspc(N,bra,ket) &
        + CGbra6*CGket6*Hspu(N,SerBra6,SerKet6)

    end do
    end do

    DIM=CBasisDim(N)
    allocate(Htmp(DIM,DIM))
    allocate(Etmp(DIM))
    allocate(Stmp(DIM,DIM))

    Htmp=0
    do bra=1, CBasisDim(N)
    do ket=1, CBasisDim(N)
        Htmp(bra,ket)=Hspc(N,bra,ket)
    end do  
    end do

    call diagonalise(DIM,DIM,EV,Htmp,Etmp,Stmp,IROT,B,Z)

    do count=1, DIM
        index=index+1
        SPE(index)=Etmp(count)
        SPBasis(index,:)=CBasis(N,count,:)
        SPWF(:,index)=Stmp(:,count)
    end do

    deallocate(Stmp)
    deallocate(Etmp)
    deallocate(Htmp)
!!! turn uncoupled Nilsson basis into coupled Nilsson basis end !!!

end subroutine NilSPE