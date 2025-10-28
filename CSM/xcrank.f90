subroutine xcrank

    use V_proce
    use V_cons
    use V_para
    use V_hom
    use V_spe
    use V_crank

    implicit none

    integer :: NF !-1 factor
    integer :: SeBra1, SeBra2, SeKet1, SeKet2
    !serial number of coupled basis in one certain cranking basis
    real :: Nb, Lb, Lamb, Sigb, Omeb
    real :: Nk, Lk, Lamk, Sigk, Omek

    integer :: DIM
    integer :: EV=1, IROT
    integer, dimension(1000,1000) :: B,Z

    real, allocatable, dimension(:) :: Etmp
    real, allocatable, dimension(:,:) :: Htmp
    real, allocatable, dimension(:,:) :: Stmp

    DIM=CrankDim(N)

    signature=-1.

    do bra=1, DIM
    do ket=bra, DIM

        SeBra1=int(CrankBV(N,bra,1))
        SeBra2=int(CrankBV(N,bra,2))
        SeKet1=int(CrankBV(N,ket,1))
        SeKet2=int(CrankBV(N,ket,2))

        Nb=UBasis(N,SeBra1,1)
        Lb=UBasis(N,SeBra1,2)
        Lamb=UBasis(N,SeBra1,3)
        Sigb=UBasis(N,SeBra1,4)
        Omeb=UBasis(N,SeBra1,5)

        Nk=UBasis(N,SeKet1,1)
        Lk=UBasis(N,SeKet1,2)
        Lamk=UBasis(N,SeKet1,3)
        Sigk=UBasis(N,SeKet1,4)
        Omek=UBasis(N,SeKet1,5)

        NF=((-1)**(Nk+0.5-signature/2))

        Hcrk(N,bra,ket)=Hspu(N,SeBra1,SeKet1) &
        + signature*((-1)**Nk)*Hspu(N,SeBra1,SeKet2) &
        + signature*((-1)**Nb)*Hspu(N,SeBra2,SeKet1) &
        + Hspu(N,SeBra2,SeKet2)

        Hcrk(N,bra,ket)=0.5*Hcrk(N,bra,ket)

        if (Lb/=Lk) cycle

        if (abs(Omeb)==abs(Omek) .and. abs(Omeb)==0.5) then
            if (Sigb==-Sigk) then
                if (Lamb==-Lamk+1) &
                    Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
                    -0.5*omega*NF*sqrt(Lk*(Lk+1)+Lamk*(-Lamk+1))
                if (Lamb==-Lamb-1) &
                    Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
                    -0.5*omega*NF*sqrt(Lk*(Lk+1)+Lamk*(-Lamk-1))
            end if

            if (Sigb==-Sigk+1) then
                if (Lamb==-Lamk) &
                    Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
                    -0.5*omega*NF*sqrt(0.75+Sigk*(-Sigk+1))
            end if

            if (Sigb==-Sigk-1) then
                if (Lamb==-Lamk) &
                    Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
                    -0.5*omega*NF*sqrt(0.75+Sigk*(-Sigk-1))
            end if
        else
            if (Sigb==Sigk) then
                if (Lamb==Lamk+1) &
                    Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
                    -0.5*omega*sqrt((Lk*(Lk+1)-(Lamk*(Lamk+1))))
                if (Lamb==Lamk-1) &
                    Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
                    -0.5*omega*sqrt((Lk*(Lk+1)-(Lamk*(Lamk-1))))
            end if

            if (Sigb==Sigk+1) then
                if (Lamb==Lamk) &
                    Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
                    -0.5*omega*sqrt(0.75-Sigk*(Sigk+1))
            end if

            if (Sigb==Sigk-1) then
                if (Lamb==Lamk) &
                    Hcrk(N,bra,ket)=Hcrk(N,bra,ket) &
                    -0.5*omega*sqrt(0.75-Sigk*(Sigk-1))
            end if
        end if

    end do
    end do

    allocate(Htmp(DIM,DIM))
    allocate(Etmp(DIM))
    allocate(Stmp(DIM,DIM))

    do bra=1, DIM
    do ket=1, DIM
        Htmp(bra,ket)=Hcrk(N,bra,ket)
    end do
    end do

    do bra=1, DIM
        write(*,'(1000f12.6)') Htmp(bra,:)
    end do

    stop

    call diagonalise(DIM,DIM,EV,Htmp,Etmp,Stmp,IROT,B,Z)

    do count=1, DIM
        write(*,*) Etmp(count)/hbom0
    end do

    do count=1, CrankDim(N)
        index=index+1
        CrkE(index)=Etmp(count)
        CrkBasis(index,:)=CrankBV(N,count,:)
        CrkWF(:,index)=Stmp(:,count)
    end do

    deallocate(Stmp)
    deallocate(Etmp)
    deallocate(Htmp)

end subroutine xcrank