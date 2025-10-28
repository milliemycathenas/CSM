module V_proce
!process variables

    integer :: count, serial, index
    integer :: bra, ket !especially for basis vectors
    integer :: braI, ketI, braII, ketII
    !especially for HFBC matrix with DIM*2 dimensions 
    integer :: DIM, ADIM !especially for matrix dimension, ADIM=dimension for all N shells
    integer :: N !especially for N shell number
    integer :: PorN !especially for two cycles of proton=1 and neutron=2
    integer :: iter !especially for interation times
    integer :: flag !flag for conditional judgement

    real :: sum !especially for add-up

end module V_proce

module V_cons
!constant or coefficient variables 
    real :: d40, d42, d44 !coefficients for hexadecapole deformation

    double precision :: pi=3.14159265359
end module V_cons

module V_para
!input parameter variables

    integer :: NL, NH !lowest and highest N shell in the valence space
    !!!N start from 0
    integer :: NuMass, ProNum, NeuNum
    !A, Z, N of the target nucleus

    real :: kappa, mu 
    !shell correlation parameter for proton and neutron

    real :: beta2, beta4, gamma, gamrad
    real :: eps2, eps22, eps4 !\epsilon_{2,22,4}
    !deformation parameters, gamrad is gamma in radian system
    
    real :: hbom0p, hbom0n !\hbar\omega_0
    real :: hbom0tp, hbom0tn, hbom0dt 
    real :: hbom0, hbom0t

    real :: GStren, GStrenP, GStrenN !pairing strength G

    real :: theta, phi, thetar, phir 
    !axis tilted angles and corresponding values in radian system
    real :: omega !rotation frequency
    real :: omx, omy, omz !\omega_x,y,z
    !omx=\omega\sin\theta\cos\phi
    !omy=\omega\sin\theta\sin\phi
    !omz=\omega\cos\theta
    
    real, dimension(10) :: kappap, mup, kappan, mun
    !shell correlation parameter for proton and neutron for each N shell
    data kappap /0.120, 0.120, 0.105, 0.090, &
    0.065, 0.060, 0.054, 0.054, 0.054, 0./
    data mup /0.00, 0.00, 0.00, 0.30, &
    0.57, 0.65, 0.69, 0.69, 0.69, 0./
    data kappan /0.120, 0.120, 0.105, 0.090, &
    0.070, 0.062, 0.062, 0.062, 0.062, 0./
    data mun /0.00, 0.00, 0.00, 0.25, &
    0.39, 0.43, 0.34, 0.26, 0.26, 0./

end module V_para

module V_dspc
!deformed cranking single-particle energy variables
    real, allocatable, dimension(:,:) :: DSPE 
    !deformed single-particle energies for each N shell
    real, allocatable, dimension(:,:) :: CrkSPE
    !cranking single-particle energies for each N shell

    real, allocatable, dimension(:,:,:) :: UBasis, CBasis
    !uncoupled basis |Nl\Lambda\Sigma> and coupled basis |Nljm>
    !UBasis(N,serial,1-4)=UBasis(N shell,serial number,1-l.2-\Lambda.3-\Sigma.4-\Omega)
    !CBasis(N,serial,1-3)=CBasis(N shell,serial number,1-l.2-j.3-m)
    !CBasis(N,serial,4-5)=index number connecting coupled and uncoupled basis
    !CBasis(N,serial,6)=1(normal state); =-1(time-reversal state)

    complex, allocatable, dimension(:,:,:) :: CrkWF
    !cranking single-particle wave functions for each N shell
end module V_dspc

module V_hfbc
!HFB cranking variables including BCS transformation
    integer :: Nval !number of valence nucleons outside of frozen core
    real :: Efermi !Fermi energy
    real :: Delta !initial pairing gap
    
    real, allocatable, dimension(:,:) :: QPE !quasi-particle energies

    complex, allocatable, dimension(:,:,:) :: QWF
    !quasi-particle wave functions for each N shell, QWF(N,basis,state)

end module V_hfbc

module V_fval
!variables to find the values of consistent angles or minima deformation
    real :: Ibeta2, Fbeta2, Ibeta4, Fbeta4
    real :: Igamma, Fgamma, Igamrad, Fgamrad
    real :: Itheta, Ftheta, Iphi, Fphi
    real :: Ithetar, Fthetar, Iphir, Fphir
    !I=initial, F=final, r=radian system


end module V_fval

module V_obs
!observable or mechanical quantity variables
    integer, allocatable, dimension(:) :: OCC
    !occupation, OCC(N,serial=1 to (N+1)(N+2))
    !for ground state, the first OCC(1-Nval)=1, OCC(Nval+1-(N+1)(N+2))=0
    !excite one particle from i-th state to k-th state, OCC(i)-1, OCC(k)+1

end module V_obs
