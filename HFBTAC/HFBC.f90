subroutine HFBC
!HFB cranking iterative calculations

    use V_proce
    use V_cons
    use V_para
    use V_dspc
    use V_hfbc

    implicit none

    real :: SPE !process variable for temporary value E'
    real :: JY !for \omega_yj_y matrix element in HFB
    real :: Lbra, Jbra, Mbra, Lket, Jket, Mket, Fbra, Fket 
    !quantum numbers for each cycle, F stands for flag (normal/time-reversal)
    real, parameter :: ATTEN=0.6 !attenuation coefficient for faster convergence

    real, allocatable, dimension(:) :: u,v !initial u_j and v_j from BCS transformation
    real, allocatable, dimension(:) :: QE, LQE !quasi-particle energies, L=last iteration

    integer, allocatable, dimension(:) :: NewIndex
    !NewIndex(i)=j means the i-th old state is the most similar to the j-th new state
    real, allocatable, dimension(:,:) :: WFoverlap !overlap of old and new U, V matrix
    !to make sure the structure of U, V matrix stays the same between each iteration

    complex, allocatable, dimension(:,:) :: MatrixU, MatrixV !U, V matrix
    complex, allocatable, dimension(:,:) :: MatrixR, MatrixK !\rho, \kappa matrix
    complex, allocatable, dimension(:,:) :: MatrixG, MatrixD, LMatrixG, LMatrixD 
    !\Gamma, \Delta matrix, L=last iteration

    complex, allocatable, dimension(:,:) :: HC !HFB cranking Hamiltonian
    complex, allocatable, dimension(:) :: EC !HFB cranking quasi-particle energies
    complex, allocatable, dimension(:,:) :: WFC !HFB quasi-particle wave functions
    complex, allocatable, dimension(:,:) :: LWFC
    !quasi-particle wave functions from last iteration
    !only wave functions corresponding to positive energies will be stored
    complex, allocatable, dimension(:,:) :: TWFC, PWFC !P=positive
    !temporary quasi-particle wave functions for overlap calculations
    !only wave functions corresponding to positive energies will be stored
    
    !!!allocate block start!!!
        allocate(u(DIM))
        allocate(v(DIM))
        allocate(QE(2*DIM))
        allocate(LQE(2*DIM))
        allocate(NewIndex(DIM))
        allocate(MatrixU(DIM,DIM))
        allocate(MatrixV(DIM,DIM))
        allocate(WFoverlap(DIM,DIM))
        allocate(MatrixG(DIM,DIM))
        allocate(MatrixD(DIM,DIM))
        allocate(LMatrixG(DIM,DIM))
        allocate(LMatrixD(DIM,DIM))
        allocate(HC(2*DIM,2*DIM))
        allocate(EC(2*DIM))
        allocate(WFC(2*DIM,2*DIM))
        allocate(LWFC(2*DIM,DIM))
        allocate(PWFC(2*DIM,DIM))
        allocate(TWFC(2*DIM,DIM))
    !!!allocate block end!!!

    if (PorN==1) GStren=GStrenP
    if (PorN==2) GStren=GStrenN

    omx=omega*sin(thetar)*cos(phir)
    omy=omega*sin(thetar)*sin(phir)
    omz=omega*cos(thetar)


    !!!initialising start!!!
        u=0
        v=0
        iter=0
        EC=0
        QE=0
        LQE=0
        WFC=0
        LWFC=0
        MatrixU=0
        MatrixV=0
        MatrixR=0
        MatrixK=0
        MatrixG=0
        MatrixD=0
    !!!initialising end!!!
    
    do count=1, DIM
        SPE=DSPE(N,count)-Efermi
        v(count)=1-SPE/sqrt(SPE**2+Delta**2)
        v(count)=0.5*v(count)
        u(count)=1-v(count)
        v(count)=sqrt(v(count))
        u(count)=sqrt(u(count))
    end do

    do bra=1, DIM
    do ket=1, DIM
        Lbra=CBasis(N,bra,1)
        Jbra=CBasis(N,bra,2)
        Mbra=CBasis(N,bra,3)
        Fbra=CBasis(N,bra,6)

        Lket=CBasis(N,ket,1)
        Jket=CBasis(N,ket,2)
        Mket=CBasis(N,ket,3)

        if (Lbra/=Lket .or. Jbra/=Jket) cycle

        if (Mbra==Mket) MatrixU(bra,ket)=u(bra)

        if (Mbra==-Mket .and. Fbra==1) MatrixV(bra,ket)=v(bra)
        if (Mbra==-Mket .and. Fbra==-1) MatrixV(bra,ket)=-v(bra) 
    end do
    end do

    101 continue
    iter=iter+1

    if (iter>100) then
        write(*,*) 'HFBC: failed to converge'
        stop
    end if

    if (iter==1) then
        MatrixG=matmul(conjg(MatrixV),transpose(MatrixV))
        MatrixD=matmul(conjg(MatrixU),transpose(MatrixV))
    end if

    LMatrixG=MatrixG
    LMatrixD=MatrixD

    MatrixG=0
    MatrixD=0
    MatrixG=ATTEN*matmul(conjg(MatrixV),transpose(MatrixV))+(1-ATTEN)*LMatrixG
    MatrixD=ATTEN*matmul(conjg(MatrixU),transpose(MatrixV))+(1-ATTEN)*LMatrixD

    HC=0
    do braI=1, DIM
        braII=braI+DIM
        HC(braI,braI)=DSPE(N,braI)-Efermi
    do ketI=1, DIM
        ketII=ketI+DIM

        HC(braI,ketI)=HC(braI,ketI)+MatrixG(braI,ketI)
        HC(braI,ketII)=MatrixD(braI,ketI)

        Lbra=CBasis(N,braI,1)
        Jbra=CBasis(N,braI,2)
        Mbra=CBasis(N,braI,3)

        Lket=CBasis(N,ketI,1)
        Jket=CBasis(N,ketI,2)
        Mket=CBasis(N,ketI,3)

        if (Lbra/=Lket .or. Jbra/=Jket) cycle
        
        JY=0
        if (Mbra==Mket+1) then
            HC(braI,ketI)=HC(braI,ketI) &
            - 0.5*sqrt(Jket*(Jket+1)-Mket*(Mket+1))*omx
            JY=JY+0.5*sqrt(Jket*(Jket+1)-Mket*(Mket+1))*omy
        else if (Mbra==Mket-1) then
            HC(braI,ketI)=HC(braI,ketI) &
            - 0.5*sqrt(Jket*(Jket+1)-Mket*(Mket-1))*omx
            JY=JY-0.5*sqrt(Jket*(Jket+1)-Mket*(Mket-1))*omy
        else if (Mbra==Mket) then
            HC(braI,ketI)=HC(braI,ketI)-Mket*omz
        end if

        HC(braI,ketI)=HC(braI,ketI)+cmplx(0.,JY)
        HC(braII,ketII)=-conjg(HC(braI,ketI))
        HC(braII,ketI)=-conjg(HC(braI,ketII))
        
    end do
    end do

    ! if (N==3) then
    ! 	do bra=1, 2*DIM
    ! 		write(*,'(1000f12.6)') real(HC(bra,1:2*DIM))
    ! 	end do
    ! end if

    call diagonalise_cmplx(2*DIM,HC,WFC,EC)

    do count=1, 2*DIM
        QE(count)=real(EC(count))
    end do

    call orderE(2*DIM, QE)

    ! write(*,*) 'iter=', iter
    ! do count=1, 2*DIM
    ! 	write(*,'(10f12.6)') LQE(count), QE(count), abs(LQE(count)-QE(count))
    ! end do

    flag=1
    do count=1, 2*DIM
        if (abs(LQE(count)-QE(count)) > 1e-4) flag=0
    end do
    if (flag==1) goto 102

    do count=1, 2*DIM
        LQE(count)=QE(count)
    end do

    if (iter==1) then
        ket=0
        do index=1, 2*DIM
            if (real(EC(index))<0) cycle
            ket=ket+1
            do braI=1, DIM
                braII=braI+DIM
                MatrixU(braI,ket)=WFC(braI,index)
                MatrixV(braI,ket)=WFC(braII,index)
            end do
            do bra=1, 2*DIM
                LWFC(bra,ket)=WFC(bra,index)
            end do
        end do		

        goto 101
    end if

    ket=0
    do count=1, 2*DIM
        if (real(EC(count))<0) cycle
        ket=ket+1
        do bra=1, 2*DIM
            PWFC(bra,ket)=WFC(bra,count)
        end do
    end do

    do bra=1, DIM
    do ket=1, DIM
        WFoverlap(bra,ket)=0
        do index=1, 2*DIM
            WFoverlap(bra,ket)=WFoverlap(bra,ket) &
            +abs(LWFC(index,bra)*PWFC(index,ket))
        end do
    end do
    end do

    call hungarian_max(DIM, WFoverlap, NewIndex)

    do ket=1, DIM
        do bra=1, 2*DIM
            TWFC(bra,ket)=PWFC(bra,NewIndex(ket))
        end do
    end do

    do ket=1, DIM
        do braI=1, DIM
            braII=braI+DIM
            MatrixU(braI,ket)=TWFC(braI,ket)
            MatrixV(braI,ket)=TWFC(braII,ket)
        end do
        do bra=1, 2*DIM
            LWFC(bra,ket)=TWFC(bra,ket)
        end do
    end do

    goto 101

    102 continue

    write(*,*) 'iteration=', iter
    do count=1, 2*DIM
        write(*,*) real(EC(count))
    end do

    do bra=1, 2*DIM
        QPE(N,bra)=real(EC(bra))
        do ket=1, DIM
            QWF(N,bra,ket)=WFC(bra,ket)
        end do
    end do

    !!!deallocate block start!!!
        deallocate(TWFC)
        deallocate(PWFC)
        deallocate(LWFC)
        deallocate(WFC)
        deallocate(EC)
        deallocate(HC)
        deallocate(LMatrixD)
        deallocate(LMatrixG)
        deallocate(MatrixD)
        deallocate(MatrixG)
        deallocate(WFoverlap)
        deallocate(MatrixV)
        deallocate(MatrixU)
        deallocate(NewIndex)
        deallocate(LQE)
        deallocate(QE)
        deallocate(v)
        deallocate(u)
    !!!deallocate block end!!!

end subroutine HFBC