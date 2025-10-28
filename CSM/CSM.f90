program CSM

    use V_proce
    use V_cons
    use V_para
    use V_hom
    use V_spe
    use V_quasi
    use V_crank

    implicit none

    integer :: DIM

    call parameters

    open(unit=101, file='omega')
    omega=-0.01
    101 omega=omega+0.01
    ! omega=0.2

    call HbarOmega

    DIM=(NH+1)*(NH+2)
    allocate(Hspu(10,DIM,DIM))
    allocate(Hspc(10,DIM,DIM))
    allocate(UBasis(10,DIM,10))
    allocate(CBasis(10,DIM,10))
    allocate(CrankBV(10,DIM,10))
    allocate(Hcrk(10,DIM,DIM))
    

    DIM=0
    do N=NL, NH
        DIM=DIM+(N+1)*(N+2)
    end do
    allocate(SPBasis(DIM,10))
    allocate(SPE(DIM))
    allocate(SPWF(DIM,DIM))
    allocate(QPE(DIM))
    allocate(quasiU(DIM))
    allocate(quasiV(DIM))
    allocate(CrkBasis(DIM/2,10))
    allocate(CrkE(DIM))
    allocate(CrkWF(DIM,DIM))
    
    call dbasis

    index=0
    Hspu=0
    Hspc=0
    do N=NL,NH
        call NilSPE
    end do
    SpeDim=index

    ! call order(SpeDim, SPE, SPBasis, SPWF)

    ! write(*,*) 'single particle energies'
    ! do count=1, SpeDim
    !     write(*,'(f12.6, 4f6.1)') SPE(count)/hbom0, SPBasis(count,1:4)
    ! end do

    ! call quasi

    ! call DpOrder(SpeDim, QPE, SPBasis, SPWF)

    ! write(*,*) 'quasi particle energies'
    ! do count=1, SpeDim
    !     write(*,'(f12.6, 4f6.1)') QPE(count)/hbom0, SPBasis(count,1:4)
    ! end do

    call crbasis

    ! index=0
    ! Hcrk=0
    ! CrkE=0
    ! CrkWF=0
    ! do N=NL,NH
    !     call QuasiCrank
    ! end do
    ! DIM=index

    ! call order(DIM, CrkE, CBasis, CrkWF)

    ! do count=1, DIM
    !     write(*,*) CrkE(count)/hbom0
    ! end do 

    ! stop

    ! write(101, '(1000f12.6)') omega/hbom0, CrkE(:)/hbom0

    ! if (abs(omega/hbom0-0.1)<=0.0001) then
    !     do count=1, DIM
    !         write(*,*) CrkE(count)/hbom0-6.485
    !     end do
    ! end if

    index=0
    Hcrk=0
    do N=NL, NH
        call xcrank
    end do
    CrkDim=index

    ! call order(CrkDim, CrkE, CrkBasis, CrkWF)

    ! write(*,*) 'cranking single particle energy'
    ! do count=CrkDim, 1, -1
    !     write(*,'(f12.4, 4f6.1)') CrkE(count)/hbom0, CrkBasis(count,1:2)
    ! end do

    stop

    ! write(101,'(1000f12.6)') omega, CrkE(:)/hbom0

    deallocate(CrkWF)
    deallocate(CrkE)
    deallocate(CrkBasis)
    deallocate(quasiV)
    deallocate(quasiU)
    deallocate(QPE)
    deallocate(SPWF)
    deallocate(SPE)
    deallocate(SPBasis)

    deallocate(Hcrk)
    deallocate(CrankBV)
    deallocate(CBasis)
    deallocate(UBasis)
    deallocate(Hspc)
    deallocate(Hspu)

    if (omega>=0.99) goto 102
    goto 101
    102 continue
    close(unit=101)

end program CSM