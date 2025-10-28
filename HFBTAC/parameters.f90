subroutine parameters

    use V_para
    use V_cons
    use V_hfbc

    implicit none

    NL=2
    NH=3

    NuMass=27
    ProNum=12
    NeuNum=NuMass-ProNum

    beta2=0.28
    beta4=0.03
    gamma=20.

    Efermi=42.
    Delta=0.75

    GStrenP=1.04
    GStrenN=0.72

    omega=0.
    theta=45
    phi=45
    thetar=theta*pi/180.
    phir=phi*pi/180.

end subroutine parameters