subroutine quasi

    use V_proce
    use V_cons
    use V_para
    use V_hom
    use V_spe
    use V_quasi

    implicit none

    if (PorN==1) then
        Nval=real(ProNum)-(NL*(NL+1)*(NL+2))/3.
    else if (PorN==0) then
        Nval=real(NeuNum)-(NL*(NL+1)*(NL+2))/3.
    end if

    if (Nval<=0) then
        write(*,*) 'QuasiTrans: the lowest N shell is too high.'
        write(*,*) 'QuasiTrans: please change NL input.'
        stop
    end if

    if (mod(Nval,2.)==0) then
        Efermi=0.5*(SPE(Nval/2)+SPE(Nval/2+1))
    else 
        Efermi=SPE(Nval/2+1)
    end if

    count=0
    do while (.true.)
        do serial=1, SpeDim
            QPE(serial)=sqrt((SPE(serial)-Efermi)**2+Delta**2)
            quasiU(serial)=0.5*(1+(SPE(serial)-Efermi)/SPE(serial))
            quasiV(serial)=0.5*(1-(SPE(serial)-Efermi)/SPE(serial))
        end do
        Nv=0
        do serial=1, SpeDim
            Nv=Nv+quasiV(serial)**2
        end do
        Nv=Nv*2
        if (abs(Nval-Nv)<1.e-8) goto 101
        Efermi=Efermi+(Nval-Nv)/Efermi

        count=count+1
        if (count>3000) then
            write(*,*) 'QuasiTrans: fermi energy calculation failed to converge.'
            stop
        end if
    end do

    101 continue

    do count=1, SpeDim
        if(SPE(count)<Efermi) QPE(count)=-QPE(count)
    end do

end subroutine quasi