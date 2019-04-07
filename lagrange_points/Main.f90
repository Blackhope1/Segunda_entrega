!  Parte2.f90 
!

!****************************************************************************
!
!  PROGRAM: Segunda entrega 
!
!  PURPOSE:  Cálculo de los puntos de Lagrange y análisis de su estabilidad (Tierra-Luna adimensional)
!
!****************************************************************************

    program Parte2

    use Diff_ope
    use Non_Linear_Systems
    use Cauchy_Problem
    use dislin_mod
    use Numerical_Recipes
    use wrapper
    
    implicit none
    
    integer,parameter :: points = 40000 !Pasos a integrar en los últimos apartados
    
    real :: x1(6), x2(6), x3(6), x4(6), x5(6), U1(6,6), U2(6,6), U3(6,6), U4(6,6), U5(6,6)
    complex :: lambda1(6), lambda2(6), lambda3(6), lambda4(6), lambda5(6), lambda6(6)
    real :: A1(6,6), A2(6,6), A3(6,6), A4(6,6), A5(6,6), t
    
    real :: tf, x(points,6), y(points,6), z(points,6), mu, dist(points), time(points), dist2(points), dist3(points)
    integer :: i, N, GL=6
    
    
    mu = 1/((5.972E24+7.349E22)/7.349E22)
    
    N=100
    
    t = 0.
    
    
    !CÁLCULO DE LOS PUNTOS CRÍTICOS DEL SISTEMA
    x1 = [0.5,0.,0.,0.,0.,0.] !Primer punto de Lagrange
    call Newton (F,x1)
    x2 = [-1.25,0.,0.,0.,0.,0.] !Segundo punto de Lagrange
    call Newton (F,x2)
    x3 = [1.,0.,0.,0.,0.,0.] !Tercer punto de Lagrange
    call Newton (F,x3)
    x4 = [-0.6,0.6,0.,0.,0.,0.] !Cuarto punto de Lagrange
    call Newton (F,x4)
    x5 = [-0.6,-0.6,0.,0.,0.,0.] !Quinto punto de Lagrange
    call Newton (F,x5)
    
    
    
    print*,"!!!!!!!!!!!!!!!!!!!!"
    print*," PUNTOS DE LAGRANGE"
    print*,"!!!!!!!!!!!!!!!!!!!!"
    print*,""
    print*,"   x1"
    print*, x1(1:3)
    print*,""
    print*,"   x2"
    print*, x2(1:3)
    print*,""
    print*,"   x3"
    print*, x3(1:3)
    print*,""
    print*,"   x4"
    print*, x4(1:3)
    print*,""
    print*,"   x5"
    print*, x5(1:3)
    
    call qplsca([x1(1),x2(1),x3(1),x4(1),x5(1)],[x1(2),x2(2),x3(2),x4(2),x5(2)],5)
    
    
    
    
    !LINEALIZACIÓN DEL SISTEMA EN LOS PUNTOS CRÍTICOS de F_linear
    call System_matrix( F_linear, x1, t, A1)
    call System_matrix( F_linear, x2, t, A2)
    call System_matrix( F_linear, x3, t, A3)
    call System_matrix( F_linear, x4, t, A4)
    call System_matrix( F_linear, x5, t, A5)    
    
    
    !CÁLCULO de AUTOVALORES DEL SISTEMA LINEALIZADO
    call Eigenvalues_QR(A1, lambda1)
    call Eigenvalues_QR(A2, lambda2)
    call Eigenvalues_QR(A3, lambda3)
    call Eigenvalues_QR(A4, lambda4)
    call Eigenvalues_QR(A5, lambda5)
    
    print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print*," AUTOVALORES DE CADA PUNTO DE LAGRANGE"
    print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    
    print*,""
    print*,"L1"
    print*, lambda1
    print*,""
    print*,"L2"
    print*, lambda2
    print*,""
    print*,"L3"
    print*, lambda3
    print*,""
    print*,"L4"
    print*, lambda4
    print*,""
    print*,"L5"
    print*, lambda5

    
    
    !PROPAGACIÓN CON PERTURBACIÓN EN L4
    
    
    
    x(1,:) = x4 + [2.6d-5,2.6d-5,0.,0.,0.,0.] !L4 + perturbación
    
    
    dist(1) = sqrt( ( x(1,1) + mu )**2 + ( x(1,2) )**2 )
    
    
    do i = 2,points
        x(i,:) = x(i-1,:)
        
        time(i) = t
        
        tf = .1*dble(i)
        
        call propagator( "DOP853", "non_linear", GL, t, x(i,:), tf) !DOP853/ODEX/ODE113
        
        t = .1*dble(i)
        
        dist(i) = sqrt( ( x(i,1) + mu )**2 + ( x(i,2) )**2 )
        
    enddo
    
    
    
    
    
    call plot(x(:,1),x(:,2), "red", "X", "Y", "", "DOPRI853, no lineal", "", .false.)
    
    call plot(time,dist, "green", "", "", "", "Distancia a la masa primaria en funcion del tiempo", "", .false.)
    
  end program Parte2

