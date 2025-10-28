subroutine HbarOmega

    use V_proce
    use V_cons
    use V_para
    
    implicit none

    integer, parameter :: dpk = kind(1.0d0)
    integer :: Loop1, Loop2, LoopX, LoopPhi
    integer, parameter :: NGauss=11
    integer :: NsecX, NsecY !number of section for subsection integral
    real :: x, DeltaX !variable substitution x=cos\theta
    real :: aphi, DeltaPhi
    real :: ITCX, ITCY !ITC=intercept coefficient (a+b)/2 from variable substitution
    real :: fx, trial, result, macroshape
    real :: LengthX, LengthY
    !range of integration of each subsection
    real(dpk), dimension(20) :: LX, LW
    

    call GaussLegendre(11, LX, LW)
    !zeros and weight for 11-point Gauss-Legendre integral

    gamrad=gamma*pi/180.
    eps2=beta2*cos(gamrad)
    eps22=beta2*sin(gamrad)/sqrt(0.75)
    eps4=beta4

    d40=(5.*cos(gamrad)**2.+1.)/6.
    d42=-sqrt(30.)*sin(2.*gamrad)/12.
    d44=sqrt(70.)*(sin(gamrad)**2.)/12.

    hbom0tp=41.*(dble(NuMass)**(-1./3.))*(1-dble(NeuNum-ProNum)/(3.*dble(NuMass)))
    hbom0tn=41.*(dble(NuMass)**(-1./3.))*(1+dble(NeuNum-ProNum)/(3.*dble(NuMass)))

    if (eps4==0) then
        hbom0dt=1./(1.+eps2/3.+eps22/2.)/(1.+eps2/3.-eps22/2.)/(1.-2.*eps2/3.)
        hbom0dt=hbom0dt**(1./3.)
        goto 101
    end if

    NsecX=2
    NsecY=3
    trial=1.

    do count=1, 6

        LengthX=1./real(NsecX)
        LengthY=pi/real(NsecY)

        result=0
        do Loop2=1, NsecY
        ITCY=(2*Loop2-1)*LengthY
        do Loop1=1, NsecX
        ITCX=-1+(2*Loop1-1)*LengthX
            do LoopPhi=1, NGauss
            aphi=LengthY*LX(LoopPhi)+ITCY
            do LoopX=1, NGauss
            x=LengthX*LX(LoopX)+ITCX
                macroshape=1.-eps2/3.*(3.*x*x-1.)+eps22/2.*(1.-x*x)*cos(2.*aphi) &
                + eps4/4.*d40*(35.*x**4.-30.*x**2+3.) &
                + 5.*eps4/sqrt(10.)*d42*(7.*x**2.-1.)*(1.-x**2.)*cos(2.*aphi) &
                + 35.*eps4/2./sqrt(70.)*d44*(1.-x**2.)*(1.-x**2.)*cos(4.*aphi)
                macroshape=1./4./pi/sqrt(1.+eps2/3.+eps22/2.) &
                /sqrt(1.+eps2/3.-eps22/2.)/sqrt(1.-2.*eps2/3.) &
                /macroshape**(1.5)
                result=result+LW(LoopX)*LW(LoopPhi)*macroshape
            end do
            end do
        end do
        end do
        hbom0dt=result*LengthX*LengthY 

        hbom0dt=hbom0dt**(1./3.)
        if (abs(trial-hbom0dt) < 1.e-5) goto 101
        NsecX=NsecX*2
        NsecY=NsecY*2
        trial=hbom0dt

        if (count>5) then
            write(*,*) 'HbarOmega: integration failed to converge.'
            stop
        end if
    end do

    101 continue
 
    hbom0p=hbom0tp*hbom0dt
    hbom0n=hbom0tn*hbom0dt

end subroutine HbarOmega