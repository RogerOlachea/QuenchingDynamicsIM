module shared
implicit none
save
integer, parameter :: L = 100, conf = 1500 !Lattice size and configuration
integer, parameter :: m = 2147483647 !2^31 - 1
end module

subroutine rand(x)
use shared
implicit none
integer, parameter :: a = 1664525, c = 374223  !multiplier, increment
integer :: x !seed
x = a*x + c
if (x < 0) then
    x = x + m
end if
end subroutine

SUBROUTINE domainwalls(State, L, D)
    IMPLICIT NONE
    INTEGER :: L
    INTEGER i, D, H
    INTEGER, DIMENSION(L) :: State
    D = 0
    DO i = 1, L
        IF(i == L) THEN
            H = State(i)*State(1)
            IF (H < 0) THEN
                D = D + 1
            END IF
        ELSE
            H = State(i)*State(i+1)
            IF (H < 0) THEN
                D = D + 1
            END IF
        END IF
    END DO
END SUBROUTINE

program QuenchingDynamics
use shared
implicit none
integer :: i !configurtions loop counter
integer :: j !lattice creation loop counter
integer :: k !Glauber loops counter
integer :: t !Montecarlo time step counter
integer :: i1 !Mesaurments writting counter
integer :: s !Glauber spin state
integer :: x,y, p!random number seeds
integer, dimension(L) :: lattice, perlat,lattice0 !Lattice size
real :: seedy!probability to flip, dice for selectin spin 
!integer, parameter :: a = 1664525, c = 1013904223  !multiplier, increment
integer :: idx
integer :: dE, DW !energy change 

!Measurments variables
real, dimension(L**2) :: DomainWall, Magnetization, Persistence, Autocorrelation

open(12, file='/home/r_o_c/Modular2/IsingModel/test0.dat')
x = 16180  !seed
y = 31637
p = 76315
lattice = 0 !Initialize variable
DomainWall = 0
Magnetization = 0
Persistence = 0
Autocorrelation = 0
!Initialize configuration
do i=1, conf
    !Create spin lattice
    do j=1, L
        call rand(x)
        if(x > real(m) / 2.0) then
            lattice(j) = 1
        else 
            lattice(j) = -1
        end if
    end do
    
    perlat = 1
    lattice0 = lattice
    !Montecarlo-time-steps-----------------------------------------------------------------
    do t=1, L**2
        !Glauber-dynamics-start----------------------------------------------------------
        do k=1, L
            dE = 0
            !random picking index between 1 an L
            call rand(y)

            seedy = real(L-1)*(real(y) / real(m)) + 1.0
        
            idx = NINT(seedy)
            s = lattice(idx) !spin state at index 'idx'
        
            !Energy change calculations for periodic boundary conditions
            IF (idx == 1) THEN
                dE = 2*(lattice(L)*s + lattice(idx + 1)*s)
            ELSE IF (idx == L) THEN
                dE = 2*(lattice(idx-1)*s + lattice(1)*s)
            ELSE
                dE = 2*(lattice(idx - 1)*s + lattice(idx + 1)*s)
            END IF

            !if energy change is less than zero, we accept the flip
            if (dE < 0) then
                lattice(idx) = -s
                perlat(idx) = 0
            !if energy change is equal to zero then we accept the flip with half probability
            else if (dE == 0) then
                call rand(p)
                if (p < real(m) / 2.0) then
                    lattice(idx) = -s
                    perlat(idx) = 0
                end if
            end if
            !Glauber-dynamics-end----------------------------------------------------------

        end do
        !Measurements-----------------------------------------------------------------------
        call domainwalls(lattice, L, DW)
        DomainWall(t) = DomainWall(t) + DW
        Magnetization(t) = Magnetization(t) + ABS(SUM(lattice))
        Persistence(t) = Persistence(t) + SUM(perlat)
        Autocorrelation(t) = Autocorrelation(t) + SUM(lattice0*lattice)
    end do
end do

do i1 = 1, L**2
    write(12,*) i1, DomainWall(i1) / (L*conf), Magnetization(i1) / (L*conf), Persistence(i1), & 
    Autocorrelation(i1) / (L*conf)
end do

end program

