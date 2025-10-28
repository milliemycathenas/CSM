!MathTools include subroutine diagonalise, OrderCmplx, CGFactor, 
!LOGFACT, COMBIN 

!!!jacobi diagonalising subroutine start!!!
        SUBROUTINE diagonalise(N1,N,EV,A,D,V,IROT,B,Z)

                    !N1：所要求解的若干哈密顿量矩阵中的最大维数
                    !N：当前求解的哈密顿量矩阵的维数
                    !EV：逻辑量，取.TRUE.表示输出本征态矢，取.FALSE.表示不输出本征态矢
                    !该逻辑量在调用JACOBI的程序中也要加以说明;LOGICAL EV,可以换成一个整型量分别取1,0代替
                    !A：二维数组，输入哈密顿量矩阵元
                    !D：一维数组，输出本征值
                    !V：二维数组，输出本征态矢
                    !IROT：子程序中的一个输出参量，由具体算法而来，表示求解过程中哈密顿量变换的次数，可不必管它。
                    !B, Z：都是二维工作数组，运算过程中要用到，在调用它们的程序（可以是主程序，也可以说是某个子程序）中也要加以说明                   
        INTEGER EV
        !DIMENSION A(N1,N1),D(N1),V(N1,N1),B(N1),Z(N1)
        DIMENSION A(N1,N1)
        DIMENSION D(N1)
        DIMENSION V(N1,N1)
        DIMENSION B(N1)
        DIMENSION Z(N1)

        !DO I=1,N
        !    DO J=1,N
        !WRITE (*,*) 'A(',J,',',I,')',A(J,I)
        !    END DO
        !END DO
        !WRITE (*,*) '子程序里的求解维数',N,N1

        IF(EV==0) GOTO 10
        DO IP=1,N
        DO IQ=1,N
            IF(IP-IQ) 50,60,50
    60	    V(IP,IP)=1.
            GOTO 70
    50	    V(IP,IQ)=0.
    70      END DO
        END DO
    10	DO IP=1,N
        D(IP)=A(IP,IP)
        B(IP)=D(IP)
        Z(IP)=0.
        END DO
        IROT=0
        DO I=1,50
        SM=0.
        NM1=N-1
        DO IP=1,NM1
            IPP1=IP+1
            DO IQ=IPP1,N
            SM=SM+ABS(A(IP,IQ))
            END DO
        END DO
        IF(SM) 110,120,110
    110	  IF(I-4) 130,140,140
    130	  TRESH=0.2*SM/(FLOAT(N)*FLOAT(N))
        GOTO 150
    140	  TRESH=0.
    150	  DO IP=1,NM1
            IPP1=IP+1
            DO IQ=IPP1,N
            G=100.*ABS(A(IP,IQ))
            IF(I.GT.4.AND.ABS(D(IP))+G.EQ.ABS(D(IP))&
            .AND.ABS(D(IQ))+G.EQ.ABS(D(IQ))) GOTO 200
            IF(ABS(A(IP,IQ)).LE.TRESH) GOTO 160
            H=D(IQ)-D(IP)
            IF(ABS(H)+G.EQ.ABS(H)) GOTO 240
            THETA=0.5*H/A(IP,IQ)
            T=1./(ABS(THETA)+SQRT(1.+THETA*THETA))
            IF(THETA.LT.0.) T=-T
            GOTO 250
    240	      T=A(IP,IQ)/H
    250	      C=1./SQRT(1.+T*T)
            S=T*C
            H=T*A(IP,IQ)
            Z(IP)=Z(IP)-H
            Z(IQ)=Z(IQ)+H
            D(IP)=D(IP)-H
            D(IQ)=D(IQ)+H
            A(IP,IQ)=0.
            IPM1=IP-1
            IF(IPM1) 260,260,270
    270	      DO J=1,IPM1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=C*G-S*H
                A(J,IQ)=S*G+C*H
            END DO
    260	      IQM1=IQ-1
            IF(IQM1-IPP1) 300,290,290
    290	      DO J=IPP1,IQM1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=C*G-S*H
                A(J,IQ)=S*G+C*H
                END DO
    300	      IQP1=IQ+1
            IF(N-IQP1) 330,320,320
    320	      DO J=IQP1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=C*G-S*H
                A(IQ,J)=S*G+C*H
            END DO
    330	      IF(EV==0) GOTO 350
            DO J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=C*G-S*H
                V(J,IQ)=S*G+C*H
                END DO
    350	      IROT=IROT+1
            GOTO 160
    200	      A(IP,IQ)=0.
    160	    END DO
        END DO
        DO IP=1,N
            B(IP)=B(IP)+Z(IP)
            D(IP)=B(IP)
            Z(IP)=0.
        END DO
        END DO
    120	RETURN
        END SUBROUTINE diagonalise
!!!jacobi diagonalising subroutine end!!!

subroutine diagonalise_cmplx(matdim,matin,matout,mateigval)
!mkl complex matrix diagonalization
    !matdim=input matrix dimension
    !matin=input complex matrix
    !matout=eigen vector (wave function)
    !mateigval=eigen value (energies)

    integer matdim,LDA,LDVL,LDVR,LWMAX,INFO,LWORK
    complex,intent(in)::matin(matdim,matdim)
    complex,intent(out)::matout(matdim,matdim)
    complex,intent(out)::mateigval(matdim)
    REAL,allocatable::RWORK(:)
    complex,allocatable::WORK(:)
    complex,allocatable::VL(:,:)
    complex,allocatable::VR(:,:)
    LDA = matdim
    LDVL = matdim
    LDVR = matdim
    LWMAX = 2 * matdim + matdim**2
    ! write(*,*)matin
    allocate(RWORK(2*matdim))
    allocate(VL(LDVL,matdim))
    allocate(VR(LDVR, matdim))
    allocate(WORK(LWMAX))
    matout = matin

    LWORK = -1
    call cgeev( 'N', 'V', matdim, matout, LDA, mateigval, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    call cgeev( 'N', 'V', matdim, matout, LDA, mateigval, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    !    STOP
    END IF
    matout = VR  !output = right eigenvectors
    return
end subroutine diagonalise_cmplx

subroutine orderE(DIM,energy)
!order only energy

    implicit none

    integer :: i,j
    integer :: DIM

    real :: energy_tmp

    real, dimension(DIM) :: energy 

    do i=1, DIM-1
    do j=1, DIM-i
        if (energy(j) > energy(j+1)) then
            energy_tmp=energy(j)
            energy(j)=energy(j+1)
            energy(j+1)=energy_tmp
        end if 
    end do
    end do


end subroutine

subroutine orderEB(DIM,energy,basis)
!order energy and basis
  ! subroutine order(DIM,energy,basis,state)
    implicit none

    integer :: i,j
    integer :: DIM

    real :: energy_tmp
    real, dimension(10) :: basis_tmp
    ! real, dimension(DIM) :: state_tmp

    real, dimension(DIM) :: energy
    real, dimension(DIM,10) :: basis
    ! real, dimension(DIM,DIM) :: state  

    do i=1, DIM-1
    do j=1, DIM-i
        if (energy(j) > energy(j+1)) then
            energy_tmp=energy(j)
            energy(j)=energy(j+1)
            energy(j+1)=energy_tmp

            basis_tmp(:)=basis(j,:)
            basis(j,:)=basis(j+1,:)
            basis(j+1,:)=basis_tmp(:)

            ! state_tmp(:)=state(:,j)
            ! state(:,j)=state(:,j+1)
            ! state(:,j+1)=state_tmp(:)
        end if 
    end do
    end do

end subroutine orderEB


subroutine GaussLegendre(ordern, nodes, weights)
    use V_proce
    use V_cons

    implicit none
    integer, parameter :: dpk = kind(1.0d0)
    !kind parameter for double-precision numbers, double-precision kine
    integer :: i !loop index for nodes, half range due to symmetry
    integer :: j !iteration counter for Newton's method
    integer :: k !loop index for Legendre polynomial calculation
    integer :: ordern !number of Gauss–Legendre integration nodes
    integer :: max_iter !maximum number of Newton iterations
    real(dpk) :: eps, p1, p2, p3, pp, z, z1
    !eps: convergence criterion for Newton iteration
    !pi: constant \pi
    !p1, p1, p3: variables for recursive calculation of Legendre polynomials.
    !pp: derivative of Legendre polynomial at node.
    !z: current approximation of node (root).
    !z1: previous approximation of node.
    real(dpk) :: nodes(20), weights(20) !generated nodes and weights

    eps = 1.0e-14_dpk
    max_iter = 100

    do i = 1, (ordern + 1)/2
        z = cos(pi * (i - 0.25_dpk) / (ordern + 0.5_dpk))

        do j = 1, max_iter
            p1 = 1.0_dpk
            p2 = 0.0_dpk

            do k = 1, ordern
                p3 = p2
                p2 = p1
                p1 = ((2.0_dpk*k - 1.0_dpk)*z*p2 - (k - 1.0_dpk)*p3)/k
            end do

        pp = ordern*(z*p1 - p2)/(z*z - 1.0_dpk)
        z1 = z
        z = z1 - p1/pp
        if (abs(z - z1) < eps) exit
        end do

        nodes(i) = -z
        nodes(ordern + 1 - i) = z
        weights(i) = 2.0_dpk/((1.0_dpk - z*z)*pp*pp)
        weights(ordern + 1 - i) = weights(i)
    end do
    
end subroutine GaussLegendre

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
