program testpes
use heh2pes
implicit none
real*8, parameter :: pi=acos(-1.0d0)
real*8 :: capr, smlr, theta, v

capr=2.97d0 !in bohr
smlr=2.0d0  !in bohr
theta=pi    !in radian

call tot_pes(capr,smlr,theta,v)

write(*,*)capr,smlr,theta,v

!Energy (MRCI+Q+LR) -0.11488649041699457 Hartree
!Energy (FCI) -0.11487377922550927 Hartree
!Energy (FCI+LR) -0.11487377922550927 Hartree
end
