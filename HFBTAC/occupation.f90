subroutine occupation
!arrange multiparticle configuration

    use V_proce
    use V_cons
    use V_para
    use V_hfbc
    use V_obs

    implicit none

    OCC=0
    do count=1, int(Nval)
        OCC(count)=1
    end do
    do count=int(Nval)+1, ADIM
        OCC(count)=0
    end do   

end subroutine occupation