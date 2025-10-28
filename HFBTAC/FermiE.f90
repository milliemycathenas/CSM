subroutine FermiE
!particle number conserved calculations for Fermi energy

	use V_proce
	use V_cons
	use V_para
	use V_dspc
	use V_hfbc

	implicit none

	real :: Nv, Ntmp !calculated nucleon number based on N=\sum_{k>0}2v_k^2
	real :: DNmu !d(v_i^2)/d(\lambda) = (\Delta^2)/(2\varepsilon'_i^3)

	real, allocatable, dimension(:) :: Etmp
	!single-particle energies for all N shells
	real, allocatable, dimension(:) :: QEtmp
	!quasi-particle energies for all N shells
	real, allocatable, dimension(:) :: v2 !v_k^2 for N=\sum_{k>0}2v_k^2

	allocate(Etmp(ADIM))
	allocate(QEtmp(ADIM))
	allocate(v2(ADIM))

	index=0
	do N=NL,NH
		do count=1, (N+1)*(N+2), 2
			index=index+1
			Etmp(index)=DSPE(N,count)
		end do
	end do

	call orderE(ADIM,Etmp)

	if (Nval<=0) then
        write(*,*) 'FermiE: the lowest N shell is too high.'
        write(*,*) 'FermiE: please change NL input.'
        stop
    end if

	if (mod(Nval,2)==0) then
        Efermi=0.5*(Etmp(Nval/2)+Etmp(Nval/2+1))
    else 
        Efermi=Etmp(Nval/2+1)
    end if

	count=0
    do while (.true.)
        do serial=1,ADIM
            QEtmp(serial)=sqrt((Etmp(serial)-Efermi)**2+Delta**2)
            v2(serial)=0.5d0*(1.d0-(Etmp(serial)-Efermi)/QEtmp(serial))
        end do
        Nv=0.d0
		DNmu=0.d0
        do serial=1, ADIM
            Nv=Nv+v2(serial)
			DNmu=DNmu+1.d0/QEtmp(serial)/QEtmp(serial)/QEtmp(serial)
        end do 
		Nv=2.d0*Nv
		DNmu=Delta*Delta*DNmu
		if (abs(Ntmp-Nv).le.1.d-8) goto 101
        Efermi=Efermi+(Nval-Nv)/DNmu
		Ntmp=Nv
        
        count=count+1
        if (count>100) then
            write(*,*) 'FermiE: fermi energy calculation failed to converge.'
            stop
        end if
    end do

    101 continue
	
	deallocate(v2)
	deallocate(QEtmp)
	deallocate(Etmp)

end subroutine FermiE