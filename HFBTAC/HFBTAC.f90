program HFBTAC

    use V_proce
    use V_cons
    use V_para
    use V_dspc
    use V_hfbc
    use V_obs

    implicit none

    integer :: DIMALL 
    !dimension of single particle energies of all N shells

    call parameters

    call HbarOmega

    DIM=(NH+1)*(NH+2)
    allocate(UBasis(10,DIM,10))
    allocate(CBasis(10,DIM,10))
    allocate(DSPE(10,DIM))
    allocate(CrkSPE(10,DIM))
    allocate(CrkWF(10,DIM,DIM))

    ADIM=0
    do N=NL,NH
        ADIM=ADIM+(N+1)*(N+2)
    end do
    allocate(OCC(ADIM))
    
    PorN=1
    if (PorN==1) Nval=real(ProNum)-(NL*(NL+1)*(NL+2))/3.
    if (PorN==0) Nval=real(NeuNum)-(NL*(NL+1)*(NL+2))/3.

    DSPE=0
    LoopDefSP: &
    do N=NL, NH

        DIM=(N+1)*(N+2)
        call basis
        call DefSP
        call SPcrank

        ! stop

    end do LoopDefSP

    call FermiE
    call occupation

    call AMSP

    DIM=(NH+1)*(NH+2)
    allocate(QPE(10,2*DIM))
    allocate(QWF(10,2*DIM,2*DIM))

    ! do N=NL, NH
    ! 	DIM=(N+1)*(N+2)
    ! 	call HFBC
    ! end do

    deallocate(QWF)
    deallocate(QPE)

    deallocate(OCC)
    deallocate(CrkWF)
    deallocate(CrkSPE)
    deallocate(DSPE)
    deallocate(CBasis)
    deallocate(UBasis)

end program HFBTAC