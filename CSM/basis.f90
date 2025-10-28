subroutine dbasis
!basis for deformed single particle energy

    use V_proce
    use V_cons
    use V_para
    use V_hom
    use V_spe
    use V_crank

    implicit none

    integer :: DIM
    integer :: LoopL, LoopLam, LoopSig, LoopJ, LoopOm

    do N=NL, NH

        UBasisDim(N)=(N+1)*(N+2)
        CBasisDim(N)=(N+1)*(N+2)

        !!!uncoupled basis start!!!
        if (mod(N,2)/=0) LDim=(N+1)/2
        if (mod(N,2)==0) LDim=(N+2)/2
        serial=1
        do LoopL=1, LDim
            LamDim=2*(N-(LoopL-1)*2)+1
        do LoopLam=1, LamDim
        do LoopSig=1, 2
            UBasis(N,serial,1)=N
            UBasis(N,serial,2)=N-(LoopL-1)*2
            UBasis(N,serial,3)=-UBasis(N,serial,2)+(LoopLam-1)
            if (LoopSig==1) UBasis(N,serial,4)=0.5
            if (LoopSig==2) UBasis(N,serial,4)=-0.5
            UBasis(N,serial,5)=UBasis(N,serial,3)+UBasis(N,serial,4)
            serial=serial+1
        end do
        end do
        end do

        !!!coupled basis start!!!
        serial=1
        do LoopL=1, LDim
            do LoopJ=1, -1, -2
                OmDim=2*(N-(LoopL-1)*2+real(LoopJ/2.))+1
                do LoopOm=1, OmDim
                    CBasis(N,serial,1)=N
                    CBasis(N,serial,2)=N-(LoopL-1)*2
                    CBasis(N,serial,3)=CBasis(N,serial,2)+LoopJ/2.
                    if (CBasis(N,serial,3)<0) cycle
                    CBasis(N,serial,4)=CBasis(N,serial,3)-(LoopOm-1)

                    do count=1, UBasisDim(N)
                        if (UBasis(N,count,2)/=CBasis(N,serial,2)) cycle
                        if (UBasis(N,count,5)==CBasis(N,serial,4)) then
                            if (UBasis(N,count,4)==0.5) CBasis(N,serial,5)=count
                            if (UBasis(N,count,4)==-0.5) CBasis(N,serial,6)=count
                        end if
                    end do
                    serial=serial+1
                end do
            end do 
        end do
        serial=serial-1
    end do   

end subroutine dbasis

subroutine crbasis
!basis for cranking single particle energy

    use V_proce
    use V_cons
    use V_para
    use V_hom
    use V_spe
    use V_crank

    implicit none

    integer :: DIM
    integer :: LoopL, LoopLam, LoopSig, LoopJ, LoopOm

    ! !!!select the normal states start!!!
    !     CBasis=0
    !     do N=NL, NH
    !         CBasisDim(N)=CBasisDim(N)/2
    !         serial=1
    !         do LoopL=1, LDim
    !             do LoopJ=1, -1, -2
    !                 OmDim=2*(N-(LoopL-1)*2+real(LoopJ/2.))+1
    !                 do LoopOm=1, OmDim, 2
    !                     CBasis(N,serial,1)=N
    !                     CBasis(N,serial,2)=N-(LoopL-1)*2
    !                     CBasis(N,serial,3)=CBasis(N,serial,2)+LoopJ/2.
    !                     if (CBasis(N,serial,3)<0) cycle
    !                     CBasis(N,serial,4)=CBasis(N,serial,3)-(LoopOm-1)
    !                     if (mod(int(2*CBasis(N,serial,3)/2+1),2)==0) &
    !                     CBasis(N,serial,4)=CBasis(N,serial,4)-1

    !                     ! do count=1, UBasisDim(N)
    !                     !     if (UBasis(N,count,2)/=CBasis(N,serial,2)) cycle
    !                     !     if (UBasis(N,count,5)==CBasis(N,serial,4)) then
    !                     !         if (UBasis(N,count,4)==0.5) CBasis(N,serial,5)=count
    !                     !         if (UBasis(N,count,4)==-0.5) CBasis(N,serial,6)=count
    !                     !     end if
    !                     ! end do
    !                     serial=serial+1
    !                 end do
    !             end do 
    !         end do
    !         serial=serial-1
    !     end do
    ! !!!select the normal states end!!!

    ! !!!recombinated vector with coupled basis!!!
    !     do N=NL, NH
    !         index=0
    !         do count=1, CBasisDim(N)
    !             if (CBasis(N,count,4)<0) cycle
    !         do serial=1, CBasisDim(N)
    !             if (CBasis(N,serial,3)==CBasis(N,count,3) &
    !             .and.CBasis(N,serial,4)==-CBasis(N,count,4)) then
    !                 index=index+1
    !                 if (mod(2*int(CBasis(N,count,4))+1,4)==1) then
    !                     CrankBV(N,index,1)=count
    !                     CrankBV(N,index,2)=serial
    !                 else
    !                     CrankBV(N,index,1)=serial
    !                     CrankBV(N,index,2)=count
    !                 end if
    !             end if
    !         end do
    !         end do
    !         CrankDim(N)=index
    !     end do  
    ! !!!recombinated vector with coupled basis!!!

    !!!recombinated vector with uncoupled basis!!!
        do N=NL,NH
            index=0
            do count=1, UBasisDim(N)
                if (UBasis(N,count,5)<0) cycle
            do serial=1, UBasisDim(N)
                if (UBasis(N,serial,2)==UBasis(N,count,2) &
                .and. UBasis(N,serial,3)==-UBasis(N,count,3) &
                .and. UBasis(N,serial,4)==-UBasis(N,count,4)) then
                    index=index+1
                    CrankBV(N,index,1)=count
                    CrankBV(N,index,2)=serial
                end if
            end do
            end do
            CrankDim(N)=index
        end do
    !!!recombinated vector with uncoupled basis!!!

end subroutine crbasis