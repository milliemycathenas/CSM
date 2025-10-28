subroutine parameters

    use V_proce
    use V_cons
    use V_para
    use V_hom
    use V_spe
    use V_crank

    implicit none

    NL=4
    NH=4

    ProNum=38
    NeuNum=46
    NuMass=ProNum+NeuNum

    PorN=1
    InteMeth=1

    beta2=0.19
    beta4=0.
    gamma=0.
    gamrad=gamma*pi/180.

    eps2=beta2*cos(gamrad)
    eps22=beta2*sin(gamrad)/sqrt(0.75)
    eps4=beta4

    d40=(5.*cos(gamrad)**2.+1.)/6.
    d42=-sqrt(30.)*sin(2.*gamrad)/12.
    d44=sqrt(70.)*(sin(gamrad)**2.)/12.

    Delta=0.12

    omega=0.01

end subroutine parameters