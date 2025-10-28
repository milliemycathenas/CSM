subroutine HbarOmega

    use V_proce
    use V_cons
    use V_para
    use V_hom
    use V_spe
    
    implicit none

    integer :: Loop1, Loop2, LoopX, LoopPhi
    integer, parameter :: NGauss=11
    integer :: NsecX, NsecY !number of section for subsection integral
    real :: x, DeltaX !variable substitution x=cos\theta
    real :: phi, DeltaPhi
    real :: ITCX, ITCY !ITC=intercept coefficient (a+b)/2 from variable substitution
    real :: fx, trial, result
    real :: LengthX, LengthY
    !range of integration of each subsection

    real(kind=8), dimension(11) :: LX = &
    (/-0.978228658146058,-0.887062599770584,-0.730152005574050, &
    -0.519096129206878,-0.269543155952687,0.000000000000000, &
    0.269543155952687,0.519096129206878,0.730152005574050, &
    0.887062599770584,0.978228658146058/)
    real(kind=8), dimension(11) :: LW = &
    (/0.05566856711616346,0.125580369447997,0.186290210927733, &
    0.233193764591847,0.262804544509907,0.272925086777901, &
    0.262804544509907,0.233193764591847,0.186290210927733, &
    0.125580369447997,0.05566856711616346/)
    !zeros and weight for 11-point Gauss-Legendre integral

    if (InteMeth==1) goto 102

    hbom0tp=41.*(dble(NuMass)**(-1./3.))*sqrt(1.-2.*(NeuNum-ProNum)/3./dble(NuMass))
    hbom0tn=41.*(dble(NuMass)**(-1./3.))*sqrt(1.+2.*(NeuNum-ProNum)/3./dble(NuMass))

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
            phi=LengthY*LX(LoopPhi)+ITCY
            do LoopX=1, NGauss
            x=LengthX*LX(LoopX)+ITCX
                macroshape=1.-eps2/3.*(3.*x*x-1.)+eps22/2.*(1.-x*x)*cos(2.*phi) &
                + eps4/4.*d40*(35.*x**4.-30.*x**2+3.) &
                + 5.*eps4/sqrt(10.)*d42*(7.*x**2.-1.)*(1.-x**2.)*cos(2.*phi) &
                + 35.*eps4/2./sqrt(70.)*d44*(1.-x**2.)*(1.-x**2.)*cos(4.*phi)
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

    102 continue

    hbom0tp=41.*(dble(NuMass)**(-1./3.))*(1-dble(NeuNum-ProNum)/(3.*dble(NuMass)))
    hbom0tn=41.*(dble(NuMass)**(-1./3.))*(1+dble(NeuNum-ProNum)/(3.*dble(NuMass)))

    DeltaX=0.001
    result=0
    do LoopX=-1000,1000
        x=LoopX*DeltaX
        fx=1.-eps2*(3.*x**2.-1.)/3.+eps4*(35.*x**4.-30.*x**2.+3.)/4.
        fx=1./(fx**1.5)
        result=result+DeltaX*fx
    end do
    hbom0dt=result/2./(1.+2.*eps2/3.)/sqrt(1.-4.*eps2/3.)
    hbom0dt=hbom0dt**(1./3.)

    101 continue
 
    hbom0p=hbom0tp*hbom0dt
    hbom0n=hbom0tn*hbom0dt

    ! write(*,*) hbom0p, hbom0n, hbom0tp, hbom0tn, hbom0dt

end subroutine HbarOmega