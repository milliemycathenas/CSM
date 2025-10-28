module V_proce
    !variables for programming
    integer :: count, serial, index
    integer :: bra, ket
    integer :: PorN !proton(=1) or neutron(=0) label
end module V_proce

module V_cons
    !variables for constants
    real, parameter :: pi=3.1415926
    real :: a20, a22, d40, d42, d44 !expansion coefficient of spherical harmonic function
end module V_cons

module V_para
    !variables for input parameters
    integer :: NuMass, ProNum, NeuNum !mass, proton and neutron number of a certain nuclide
    integer :: N, NL, NH !shell N, lowest shell NL, highest shell NH

    real :: beta2, beta4, gamma, gamrad !deformation parameters
    real :: eps2, eps22, eps4 !eps2=beta2\cos\gamma, eps22=2*beta2\sin\gamma/\sqrt3
    real :: kappa, mu 
    !shell correlation parameter for proton and neutron
    real :: Delta !pairing gap
    real, dimension(10) :: kappap, mup, kappan, mun
    !shell correlation parameter for proton and neutron for each N shell

    complex :: imunit=(0.,1.) !imaginary unit i

    data kappap /0.120, 0.120, 0.105, 0.090, &
    0.065, 0.060, 0.054, 0.054, 0.054, 0./
    data mup /0.00, 0.00, 0.00, 0.30, &
    0.57, 0.65, 0.69, 0.69, 0.69, 0./
    data kappan /0.120, 0.120, 0.105, 0.090, &
    0.070, 0.062, 0.062, 0.062, 0.062, 0./
    data mun /0.00, 0.00, 0.00, 0.25, &
    0.39, 0.43, 0.34, 0.26, 0.26, 0./
end module V_para

module V_hom
    !variables for \hbar\omega
    integer :: InteMeth !integral method
    !InteMeth=1:rectangular; InteMeth=2:Gauss-Legendre

    real :: hbom0dt, hbom0p, hbom0n, hbom0, hbom0tp, hbom0tn, hbom0t
    !\hbar\omega_0 and \tilde{\hbar\omega_0} for proton and neutron
    real :: macroshape !the macroshape function of nuclear surface
end module V_hom

module V_spe
    !variables for single particle energy
    integer :: TUBasisDim, TCBasisdim 
    !dimension of uncoupled and coupled basis for all N shells
    integer :: LDim, LamDim, OmDim
    !number of l, \Lambda and \Omega that can be taken under a certain N 
    integer :: SpeDim !total number of single-particle energies for all N shells

    integer, dimension(10) :: UBasisDim 
    !dimension of uncoupled basis |Nl\Lambda\Sigma(\Omega)> for all N shells
    integer, dimension(10) :: CBasisDim 
    !dimension of coupled basis |Nlj\Omega> for all N shells
    real(kind=4), allocatable, dimension(:) :: Espu !single particle energy with UBasis
    real(kind=4), allocatable, dimension(:) :: Espc !single particle energy with CBasis
    real, allocatable, dimension(:) :: SPE !single particle energies for all N shells
    
    real, allocatable, dimension(:,:) :: CSPState 
    !eigen state of single particle energy under coupled basis
    real, allocatable, dimension(:,:) :: USPState 
    !eigen state of single particle energy under coupled basis
    real, allocatable, dimension(:,:) :: SPBasis
    !coupled basis for all N shells   
    real, allocatable, dimension(:,:) :: SPWF !single particle wave functions for all N shells

    real, allocatable, dimension(:,:,:) :: UBasis !uncoupled basis
    !1-N, 2-l, 3-\Lambda, 4-\Sigma, 5-\Omega
    real, allocatable, dimension(:,:,:) :: CBasis !coupled basis
    !1-N, 2-l, 3-j, 4-\Omega, 5-CG+, 6-CG-
    real(kind=4), allocatable, dimension(:,:,:) :: Hspu !single particle Hamiltonian with UBasis
    real(kind=4), allocatable, dimension(:,:,:) :: Hspc
    !single particle Hamiltonian with CBasis for all N shells
end module V_spe

module V_quasi
    !variables for quasi-particle transformation
    double precision :: Nv !N=2\sum_{\mu}V_{\mu}^2
    double precision :: Nval !valence particle number in the non-freezed N shell(s)
    double precision :: Efermi !calculated Fermi energy

    double precision, allocatable, dimension(:) :: QPE !quasi particle energy
    double precision, allocatable, dimension(:) :: quasiU, quasiV 
    !U_{\mu}^2 = \frac{1}{2}[1+\frac{e_{\mu}-\lambda}{\sqrt(e_{\mu}-\lambda)^2+\Delta^2}]
    !V_{\mu}^2 = \frac{1}{2}[1-\frac{e_{\mu}-\lambda}{\sqrt(e_{\mu}-\lambda)^2+\Delta^2}]
end module V_quasi

module V_crank
    !variables for x axis rotation
    integer :: CrkDim
    real :: omega !cranking frequency
    real :: signature !\alpha=0.5 or -0.5

    integer, dimension(10) :: CrankDim
    !dimension of cranking recombination basis
    
    real, allocatable, dimension(:) :: CrkE
    real, allocatable, dimension(:,:) :: CrkBasis
    real, allocatable, dimension(:,:) :: CrkWF

    real, allocatable, dimension(:,:,:) :: CrankBV
    !cranking recombination basis |\chi, \alpha>=[|\chi>+2\alpha(-1)^{N_{\chi}}|-\chi>]
    !CrankBV(N,count,1)=1->|Nljm>, CrankBV(N,count,2)->R_1(\pi)|Nljm>
    real, allocatable, dimension(:,:,:) :: Hcrk
    !cranking single particle hamiltonian for all shells
    
end module V_crank
