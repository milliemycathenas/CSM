subroutine DefSP

    use V_proce
    use V_cons
    use V_para
    use V_dspc

    implicit none

    !!!variables for CG factor subroutine start!!!
        double precision :: J1, J2, J3, M1, M2, M3
        double precision :: CG200,CG2mm,CG400,CG4mm
        
        double precision :: J1bra, J3bra, M3bra, M4bra, M5bra
        double precision :: J1ket, J3ket, M3ket, M4ket, M5ket
        double precision :: CGbra4, CGbra5, CGket4, CGket5
    !!!variables for CG factor subroutine end!!!

    real :: Nb,Lb,Lamb,Sigb, Nk,Lk,Lamk,Sigk
    real :: sqll !(-2/3)*\sqrt{\frac{2l+1}{2l'+1}}
    real :: NExbra, NExket !eliminate basis term with serial number 0

    real, allocatable, dimension(:) :: Etmp !eigen energy
    real, allocatable, dimension(:,:) :: Htmp, Htmpc !hamiltionian
    real, allocatable, dimension(:,:) :: Stmp !eigen state (wave function)

    integer :: EV=1, IROT
    integer, dimension(1000,1000) :: B,Z

    allocate(Etmp(DIM))
    allocate(Htmp(DIM,DIM))
    allocate(Htmpc(DIM,DIM))
    allocate(Stmp(DIM,DIM))

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

    do bra=1, DIM
    do ket=bra, DIM

        Nb=real(N)
        Lb=UBasis(N,bra,1)
        Lamb=UBasis(N,bra,2)
        Sigb=UBasis(N,bra,3)

        Nk=real(N)
        Lk=UBasis(N,ket,1)
        Lamk=UBasis(N,ket,2)
        Sigk=UBasis(N,ket,3)

        if (Lb/=Lk .and. Lb/=Lk+2 .and. Lb/=Lk-2) cycle

        J2=dble(UBasis(N,ket,1))
        J3=dble(UBasis(N,bra,1))
        M2=dble(UBasis(N,ket,2))
        M3=dble(UBasis(N,bra,2))

        sqll=sqrt((2*Lk+1)/(2*Lb+1))

        if (Sigb==Sigk) then
            call CGFactor(dble(2),J2,J3,dble(0),dble(0),dble(0), CG200)
            call CGFactor(dble(4),J2,J3,dble(0),dble(0),dble(0), CG400)
            if (Lamb==Lamk) then
                call CGFactor(dble(2),J2,J3,dble(0),M2,M3, CG2mm)
                call CGFactor(dble(4),J2,J3,dble(0),M2,M3, CG4mm)
                Htmp(bra,ket)=hbom0*sqll &
                *((-2./3.)*eps2*cos(gamrad)*CG200*CG2mm+eps4*d40*CG400*CG4mm)
            end if
            if (Lamb==Lamk+2) then
                call CGFactor(dble(2),J2,J3,dble(2),M2,M3, CG2mm)
                call CGFactor(dble(4),J2,J3,dble(2),M2,M3, CG4mm)
                Htmp(bra,ket)=hbom0*sqll &
                *((2./3.)*eps2*sin(gamrad)*CG200*CG2mm/sqrt(2.)+eps4*d42*CG400*CG4mm)
            end if
            if (Lamb==Lamk-2) then
                call CGFactor(dble(2),J2,J3,dble(-2),M2,M3, CG2mm)
                call CGFactor(dble(4),J2,J3,dble(-2),M2,M3, CG4mm)
                Htmp(bra,ket)=hbom0*sqll &
                *((2./3.)*eps2*sin(gamrad)*CG200*CG2mm/sqrt(2.)+eps4*d42*CG400*CG4mm)
            end if
            if (Lamb==Lamk+4) then
                call CGFactor(dble(4),J2,J3,dble(4),M2,M3, CG4mm)
                Htmp(bra,ket)=hbom0*sqll*eps4*d44*CG400*CG4mm
            end if
            if (Lamb==Lamk-4) then
                call CGFactor(dble(4),J2,J3,dble(-4),M2,M3, CG4mm)
                Htmp(bra,ket)=hbom0*sqll*eps2*d44*CG400*CG4mm
            end if

            if (Lb==Lk) then
                Htmp(bra,ket)=Htmp(bra,ket)*(Nk+1.5)
            else if (Lb==Lk-2) then
                Htmp(bra,ket)=Htmp(bra,ket)*sqrt(Nk-Lk+2.)*sqrt(Nk+Lk+1.)
            else if (Lb==Lk+2) then
                Htmp(bra,ket)=Htmp(bra,ket)*sqrt(Nk-Lk)*sqrt(Nk+Lk+3.)
            end if
        end if

        if (Sigb==Sigk .and. Lamb==Lamk .and. Lb==Lk) then
            Htmp(bra,ket)=Htmp(bra,ket)+hbom0*(Nk+1.5) &
            - 2.*kappa*hbom0t*Lamk*Sigk &
            - kappa*mu*hbom0t &
            *(Lk*(Lk+1.)-0.5*Nk*(Nk+3.))
        end if

        if (Sigb==Sigk+1 .and. Lamb==Lamk-1 .and. Lb==Lk) then
            Htmp(bra,ket)=-kappa*hbom0t &
            *sqrt(Lk*(Lk+1)-Lamk*(Lamk-1)) &
            *sqrt(0.5*1.5-Sigk*(Sigk+1))
        end if

        if (Sigb==Sigk-1 .and. Lamb==Lamk+1 .and. Lb==Lk) then
            Htmp(bra,ket)=-kappa*hbom0t &
            *sqrt(Lk*(Lk+1)-Lamk*(Lamk+1)) &
            *sqrt(0.5*1.5-Sigk*(Sigk-1))
        end if

        Htmp(ket,bra)=Htmp(bra,ket)

    end do
    end do

!!! turn uncoupled Nilsson basis into coupled Nilsson basis start !!!
    Htmpc=0
    do bra=1, DIM
    do ket=1, DIM

        J1bra=dble(CBasis(N,bra,1))
        J1ket=dble(CBasis(N,ket,1))
        M4bra=dble(UBasis(N,int(CBasis(N,bra,4)),2))
        M4ket=dble(UBasis(N,int(CBasis(N,ket,4)),2))
        M5bra=dble(UBasis(N,int(CBasis(N,bra,5)),2))
        M5ket=dble(UBasis(N,int(CBasis(N,ket,5)),2))
        J3bra=dble(CBasis(N,bra,2))
        J3ket=dble(CBasis(N,ket,2))
        M3bra=dble(CBasis(N,bra,3))
        M3ket=dble(CBasis(N,ket,3))

        call CGFactor(J1bra,dble(0.5),J3bra,M4bra,dble(0.5),M3bra, CGbra4)
        call CGFactor(J1bra,dble(0.5),J3bra,M5bra,dble(-0.5),M3bra, CGbra5)
        call CGFactor(J1ket,dble(0.5),J3ket,M4ket,dble(0.5),M3ket, CGket4)
        call CGFactor(J1ket,dble(0.5),J3ket,M5ket,dble(-0.5),M3ket, CGket5)

        NExbra=1
        NExket=1
        if (CBasis(N,bra,5)==0) NExbra=0
        if (CBasis(N,ket,5)==0) NExket=0

        Htmpc(bra,ket)=CGbra4*CGket4*Htmp(CBasis(N,bra,4),CBasis(N,ket,4)) &
        + CGbra4*CGket5*Htmp(CBasis(N,bra,4),CBasis(N,ket,5))*NExbra &
        + CGbra5*CGket4*Htmp(CBasis(N,bra,5),CBasis(N,ket,4)) &
        + CGbra5*CGket5*Htmp(CBasis(N,bra,5),CBasis(N,ket,5))*NExket

    end do
    end do
!!! turn uncoupled Nilsson basis into coupled Nilsson basis end !!!

    call diagonalise(DIM,DIM,EV,Htmpc,Etmp,Stmp,IROT,B,Z)

    do count=1, DIM
        DSPE(N,count)=Etmp(count)
    end do

    deallocate(Stmp)
    deallocate(Htmpc)
    deallocate(Htmp)
    deallocate(Etmp)
    
end subroutine DefSP