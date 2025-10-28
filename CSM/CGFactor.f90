subroutine CGFactor(J1,J2,J3,M1,M2,M3, CGResult)

	implicit none

	integer :: mu
	integer :: J123, J231, J312, J123p1
	!J123=j1+j2-j3, J231=j2+j3-j1, J312=j3+j1-j2, J123p1=j1+j2+j3+1
	integer :: J123mu, J1mu, J2mu, J32mu, J31mu
	!J123mu=j1+j2-j3-mu, J1mu=j1-m1-mu, J2mu=j2+m2-mu
	!J32=j3-j2+m1-mu, J31=j3-j1-m2+mu
	double precision :: LOGFACT !factorial result
	double precision :: SUM !sum for mu
	double precision :: SumFactor
	double precision :: J1, J2, J3, M1, M2, M3, CGResult

	CGResult=0
	if (M1+M2 /= M3) goto 101

	J123=int(J1+J2-J3)
	J231=int(J2+J3-J1)
	J312=int(J3+J1-J2)
	J123p1=int(J1+J2+J3+1.)

	CGResult=(2.*J3+1.)*dexp(LOGFACT(J123))*dexp(LOGFACT(J231))*dexp(LOGFACT(J312)) &
	*dexp(LOGFACT(int(J1-M1)))*dexp(LOGFACT(int(J1+M1))) &
	*dexp(LOGFACT(int(J2-M2)))*dexp(LOGFACT(int(J2+M2))) &
	*dexp(LOGFACT(int(J3-M3)))*dexp(LOGFACT(int(J3+M3))) &
	/dexp(LOGFACT(J123p1))

	CGResult=sqrt(CGResult)

	SUM=0
	do mu=0, 100
		if (J1+J2-J3-mu<0) goto 102
		if (J1-M1-mu<0) goto 102
		if (J2+M2-mu<0) goto 102
		if (J3-J2+M1+mu<0) cycle
		if (J3-J1-M2+mu<0) cycle
		J123mu=J1+J2-J3-mu
		J1mu=J1-M1-mu
		J2mu=J2+M2-mu
		J32mu=J3-J2+M1+mu
		J31mu=J3-J1-M2+mu
		SumFactor=dexp(LOGFACT(mu))*dexp(LOGFACT(J123mu))*dexp(LOGFACT(J1mu)) &
		*dexp(LOGFACT(J2mu))*dexp(LOGFACT(J32mu))*dexp(LOGFACT(J31mu))
		SumFactor=1/SumFactor
		SumFactor=SumFactor*(-1)**mu
		SUM=SUM+SumFactor
	end do

	102 continue

	CGResult=CGResult*SUM

	101 return

end subroutine CGFactor


function LOGFACT(n)
	!LOGFACT(n) calculates dlog(n!)

	implicit none
	integer :: count
	integer :: n
	double precision :: LOGFACT

	LOGFACT=dlog(dble(1))
	if (n==0 .or. n==1) goto 101

	do count=1, n
		LOGFACT=LOGFACT+dlog(dble(count))
	end do

	101 return

end function LOGFACT

function COMBIN(n,m)
	!COMBIN calculates the combination C_n^m

	implicit none
	integer :: m,n
	integer :: COMBIN
	double precision :: result
	double precision :: LOGFACT

	result=dexp(LOGFACT(n))/dexp(LOGFACT(m))/dexp(LOGFACT(n-m))

	COMBIN=int(result)+1

	return

end function COMBIN