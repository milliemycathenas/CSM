subroutine SPcrank
!single-particle cranking energy calculations

	use V_proce
	use V_cons
	use V_para
	use V_dspc
	
	implicit none

	real :: JY !for \omega_yj_y matrix element in HFB
	real :: Lbra, Jbra, Mbra, Lket, Jket, Mket
	!quantum numbers for each cycle

	complex, allocatable, dimension(:,:) :: HC !HFB cranking Hamiltonian
	complex, allocatable, dimension(:) :: EC !HFB cranking quasi-particle energies
	complex, allocatable, dimension(:,:) :: WFC !HFB quasi-particle wave functions

	allocate(HC(DIM,DIM))
	allocate(EC(DIM))
	allocate(WFC(DIM,DIM))

	omx=omega*sin(thetar)*cos(phir)
	omy=omega*sin(thetar)*sin(phir)
	omz=omega*cos(thetar)

	HC=0
	do bra=1, DIM
		HC(bra,bra)=DSPE(N,bra)
	do ket=1, DIM
		Lbra=CBasis(N,bra,1)
		Jbra=CBasis(N,bra,2)
		Mbra=CBasis(N,bra,3)

		Lket=CBasis(N,ket,1)
		Jket=CBasis(N,ket,2)
		Mket=CBasis(N,ket,3)

		if (Lbra/=Lket .or. Jbra/=Jket) cycle
		
		JY=0
		if (Mbra==Mket+1) then
			HC(bra,ket)=HC(bra,ket) &
			- 0.5*sqrt(Jket*(Jket+1)-Mket*(Mket+1))*omx
			JY=JY+0.5*sqrt(Jket*(Jket+1)-Mket*(Mket+1))*omy
		else if (Mbra==Mket-1) then
			HC(bra,ket)=HC(bra,ket) &
			- 0.5*sqrt(Jket*(Jket+1)-Mket*(Mket-1))*omx
			JY=JY-0.5*sqrt(Jket*(Jket+1)-Mket*(Mket-1))*omy
		else if (Mbra==Mket) then
			HC(bra,ket)=HC(bra,ket)-Mket*omz
		end if

		HC(bra,ket)=HC(bra,ket)+cmplx(0.,JY)
	end do
	end do

	call diagonalise_cmplx(DIM,HC,WFC,EC)

	do count=1, DIM
		CrkSPE(N,count)=real(EC(count))
	end do

	do bra=1, DIM
	do ket=1, DIM
		CrkWF(N,bra,ket)=WFC(bra,ket)
	end do
	end do

	deallocate(WFC)
	deallocate(EC)
	deallocate(HC)

end subroutine SPcrank