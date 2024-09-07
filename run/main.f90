    program testpes
    use heh2pes
    implicit none
    integer,parameter:: rk = kind(1.0d0)
    integer,parameter:: nR = 1000
    real*8, parameter :: pi=acos(-1.0d0)
    real(rk),parameter:: RMin = 1.5_rk,& 
    Rmax = 25.0_rk,&
    length = Rmax - Rmin,& 
    dR = length/real((nR+1),rk) 
    real*8 :: capr, smlr, theta
    real(rk),allocatable::Rs(:),V(:)
    integer:: i 
    ! capr=2.97d0 !in bohr
    smlr=1.4_rk!2.0d0  !in bohr
    theta=pi    !in radian
    allocate(Rs(nR),source=0.0_rk)
    allocate(V(nR),source=0.0_rk)
    do i=1,nR 
        Rs(i) = RMin+i*dR
        call tot_pes(Rs(i),smlr,theta,v(i))
    end do 

    V = V - V(nR)
    print*,'minval=',minval(v),'minloc',minloc(v),'mincoord',Rs(minloc(v))
    open(66,file='../pot_He_vs_H_2.dat')
        do i=1,nR
            write(66,'(2(g0,x))')Rs(i),V(i)
        end do 
    close(66)
    

    !Energy (MRCI+Q+LR) -0.11488649041699457 Hartree
    !Energy (FCI) -0.11487377922550927 Hartree
    !Energy (FCI+LR) -0.11487377922550927 Hartree
    end
