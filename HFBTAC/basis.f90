subroutine basis

    use V_proce
    use V_cons
    use V_para
    use V_dspc

    implicit none

    integer :: LDim, LamDim, OmDim
    !number of l, \Lambda and \Omega that can be taken under a certain N
    integer :: LoopL, LoopLam, LoopSig, LoopOm, LoopJ
    !counter for l, \Lambda, \Sigma, \Omega and J
        
    if (mod(N,2)/=0) LDim=(N+1)/2
    if (mod(N,2)==0) LDim=(N+2)/2

    serial=1
    do LoopL=1, LDim
        LamDim=2*(N-(LoopL-1)*2)+1
    do LoopLam=1, LamDim
    do LoopSig=1, 2
        UBasis(N,serial,1)=N-(LoopL-1)*2
        UBasis(N,serial,2)=-UBasis(N,serial,1)+(LoopLam-1)
        if (LoopSig==1) UBasis(N,serial,3)=0.5
        if (LoopSig==2) UBasis(N,serial,3)=-0.5
        UBasis(N,serial,4)=UBasis(N,serial,2)+UBasis(N,serial,3)
        serial=serial+1
    end do
    end do
    end do

    serial=1
    do LoopL=1, LDim
    do LoopJ=1, -1, -2
        OmDim=N-(LoopL-1)*2+LoopJ/2.+0.5
        if (mod(OmDim,2)==1) flag=1
        if (mod(OmDim,2)==0) flag=-1
    do LoopOm=1, OmDim
        CBasis(N,serial,1)=N-(LoopL-1)*2
        CBasis(N,serial,2)=CBasis(N,serial,1)+LoopJ/2.
        if (CBasis(N,serial,2)<0) cycle
        CBasis(N,serial,3)=CBasis(N,serial,2)-(LoopOm-1)*2
        do count=1, DIM
            if (UBasis(N,count,1)/=CBasis(N,serial,1)) cycle
            if (UBasis(N,count,4)==CBasis(N,serial,3)) then
                if (UBasis(N,count,3)==0.5) CBasis(N,serial,4)=count
                if (UBasis(N,count,3)==-0.5) CBasis(N,serial,5)=count
            end if
        end do
        CBasis(N,serial,6)=flag
        serial=serial+1
        CBasis(N,serial,1)=CBasis(N,serial-1,1)
        CBasis(N,serial,2)=CBasis(N,serial-1,2)
        CBasis(N,serial,3)=-CBasis(N,serial-1,3)
        do count=1, DIM
            if (UBasis(N,count,1)/=CBasis(N,serial,1)) cycle
            if (UBasis(N,count,4)==CBasis(N,serial,3)) then
                if (UBasis(N,count,3)==0.5) CBasis(N,serial,4)=count
                if (UBasis(N,count,3)==-0.5) CBasis(N,serial,5)=count
            end if
        end do
        CBasis(N,serial,6)=-flag
        serial=serial+1
    end do
    end do 
    end do
    
end subroutine basis