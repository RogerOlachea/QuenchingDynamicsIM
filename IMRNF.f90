PROGRAM test2
    IMPLICIT NONE
    INTEGER, parameter :: L = 500, conf = 1500
    REAL, PARAMETER :: pi = 3.1416
    INTEGER, DIMENSION(L) :: initState, PerLat, init
    REAL, DIMENSION(L**2) :: DW, Mag, PerVec, AC
    INTEGER :: i, D, t, k, s, j, p, idx,g
    REAL :: dE, h, X1, R, THETA

    !Random number variables
    INTEGER :: a = 1664525, c = 1013904223, m = 2147483647, x
    OPEN(12, FILE = '/home/r_o_c/Modular2/IsingModel/Data/IMNRF/IM500Data.txt')
    j = 65539
    x = 1664525
    P = 65539
    idx = 0
    h = 0

    DW = 0
    Mag = 0
    PerVec = 0
    AC = 0.0
    
    DO k = 1, conf
        !INITIAL CONFIGURATION
        DO i = 1, L
            x = a*x + c
            IF (x < 0) THEN
                x = x + m
            END IF
            
            !50% are -1, 50% are 1
            IF(x > FLOAT(m) / 2) THEN
                initState(i) = -1
            ELSE
                initState(i) = 1
            END IF

        END DO
        init = initState
        PerLat = 1
        !MONTECARLO STEPS
        DO t = 1, L**2
            !GLAUBER DYNAMICS
            DO i = 1, L-1
                dE = 0
            !Random choice of spin
                j = a*j + c
                IF (j < 0) THEN
                    j = j + m
                END IF
                
                idx = NINT((L-1)*FLOAT(j)/m + 1)
            !Random field
                g = a*g + c 
                IF (g < 0) THEN
                    g = g + m
                END IF

                X1 = (FLOAT(g) / FLOAT(m))
        
                R = SQRT(-2*LOG(X1))
                THETA = 2*pi*(99*X1 + 1)
        
                h = R*COS(THETA)
            !Spin
                s = initState(idx)
            !Energy change of flipping spin
                IF (j == 1) THEN 
                    dE = 2.0*(FLOAT(initState(L)*s) + FLOAT(initState(idx + 1)*s)) + 2.0*FLOAT(s)*h
                ELSE IF (idx == L) THEN
                    dE = 2.0*(FLOAT(initState(idx-1)*s) + FLOAT(initState(1)*s)) + 2.0*FLOAT(s)*h
                ELSE
                    dE = 2.0*(FLOAT(initState(idx - 1)*s) + FLOAT(initState(idx + 1)*s)) + 2.0*FLOAT(s)*h
                END IF
                !If energy change is less than zero, we accept the flip
                IF (dE < 0) THEN
                    initState(idx) = -s
                    PerLat(idx) = 0
                !if energy change is equal to zero then we accept the flip with half probability
                ELSE IF (dE == 0) THEN
                    p = a*p + c 
                    IF (p < 0) THEN
                        p = p + m
                    END IF
                    
                    IF (p < FLOAT(m) / 2) THEN
                        initState(idx) = -s
                        PerLat(idx) = 0
                    END IF
                END IF
            END DO
            !MEASUREMENTS

            !Domain walls
            CALL DomainWalls(initState, L, D)
            DW(t) = DW(t) + D 

            !Magnetization
            Mag(t) = Mag(t) + ABS(SUM(initState))

            !Persistance
            PerVec(t) = PerVec(t) + SUM(PerLat)

            !Autocorrelation
            AC(t) = AC(t) + SUM(init*initState)
            

        END DO
        write(*, '(A)') achar(27)//'[2K' // achar(13)
        WRITE(*,'(A, I4)', advance='no') 'MC steps: ', k
    END DO
    DO i = 1, L**2
        WRITE(12,*) i, DW(i) / (conf*L) , Mag(i)/ (conf*L), PerVec(i)/ (conf*L) , AC(i) / (conf*L)
    END DO
END PROGRAM test2

SUBROUTINE DomainWalls(State, L, D)
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
