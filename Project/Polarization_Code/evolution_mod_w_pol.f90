module evolution_mod_w_pol
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-8
  real(dp),                private :: x_init
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter, private :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:,:) :: Theta_p
  real(dp), allocatable, dimension(:,:,:) :: Nu
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta
  real(dp), allocatable, dimension(:,:,:) :: dTheta_p
  real(dp), allocatable, dimension(:,:,:) :: dNu

  ! Fourier mode list
  real(dp), allocatable, dimension(:)     :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 8+3*lmax_int

  real(dp), allocatable, dimension(:)     :: dydx


contains

  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

  end subroutine get_hires_source_function


  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i, k
    real(dp)     :: f_v

    f_v    = 0.405d0

    x_init = log(a_init)

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do k=1,n_k
       ks(k) = k_min + (k_max-k_min)*((k-1.d0)/(n_k-1.d0))**2
    end do

    ! Allocate arrays for perturbation quantities
    allocate(Theta(1:n_t, 0:lmax_int, n_k))
    allocate(Theta_p(1:n_t, 0:lmax_int, n_k))
    allocate(Nu(1:n_t, 0:lmax_int, n_k))
    allocate(delta(1:n_t, n_k))
    allocate(delta_b(1:n_t, n_k))
    allocate(v(1:n_t, n_k))
    allocate(v_b(1:n_t, n_k))
    allocate(Phi(1:n_t, n_k))
    allocate(Psi(1:n_t, n_k))
    allocate(dPhi(1:n_t, n_k))
    allocate(dPsi(1:n_t, n_k))
    allocate(dv_b(1:n_t, n_k))
    allocate(dTheta(1:n_t, 0:lmax_int, n_k))
    allocate(dTheta_p(1:n_t,0:lmax_int,n_k))
    allocate(dNu(1:n_t, 0:lmax_int, n_k))

    Theta(:,:,:)   = 0.d0
    Theta_p(:,:,:) = 0.d0
    Nu(:,:,:)       = 0.d0
    dTheta(:,:,:)  = 0.d0
    dPhi(:,:)      = 0.d0
    dPsi(:,:)      = 0.d0
    dTheta_p(:,:,:)= 0.d0
    dNu(:,:,:)      = 0.d0

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(1,:)     = 1.d0
    Psi(1,:)     = -Phi(1,:)
    delta(1,:)   = 1.5d0*Phi(1,:)
    delta_b(1,:) = delta(1,:)
    Theta(1,0,:) = 0.5d0*Phi(1,:)
       
    do i = 1, n_k
       v(1,i)         = c*ks(i)/(2*get_H_p(x_init))*Phi(1,i)
       v_b(1,i)       = v(1,i)
       Theta(1,1,i)   = -c*ks(i)/(6*get_H_p(x_init))*Phi(1,i)
       Theta(1,2,i)   = -8.d0*c*ks(i)/(15.d0*get_H_p(x_init)*get_dtau(x_init))*Theta(1,1,i)
       Theta_p(1,1,i) = -c*ks(i)/(4.d0*get_H_p(x_init)*get_dtau(x_init))*Theta(1,2,i) 
       do l = 3, lmax_int
          Theta(1,l,i) = -l/(2.d0*l+1.d0)*c*ks(i)/(get_H_p(x_init)*get_dtau(x_init))*Theta(1,l-1,i)
       end do
    end do
    
    Theta_p(1,0,:) = 1.25d0*Theta(1,2,:)
    Theta_p(1,2,:) = 0.25d0*Theta(1,2,:)

    Nu(1,0,:)       = 0.5d0*Phi(1,:)
    Nu(1,1,:)       = Theta(1,1,:)

    do i = 1,n_k
       Nu(1,2,i)    = -(c*ks(i)*a_init)**2.d0*Phi(1,i)/((12.d0*H_0**2.d0*Omega_nu)*(5.d0/(2.d0*f_v)+1))
       do l = 3, lmax_int
          Nu(1,l,i) = c*ks(i)/((2.d0*l+1)*get_H_p(x_init))*Nu(1,l-1,i)
       end do
    end do

    write(*,*) "Pertubation Equations Initialized."

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l,j_tc
    real(dp)     :: x1, x2,H,ck,ckH,a,dtau,x,bleta
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2,Pi_c

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling
    real(dp), allocatable, dimension(:) :: x_temp
    real(dp), allocatable, dimension(:) :: x_post
    real(dp), allocatable, dimension(:) :: x_total
    real(dp), allocatable, dimension(:) :: prints

    eps    = 1.d-8
    hmin   = 0.d0
    h1     = 1.d-5

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(npar))
    allocate(prints(6))

    prints(1)        = 1
    prints(2)        = 10
    prints(3)        = 30
    prints(4)        = 50
    prints(5)        = 80
    prints(6)        = 100

    y_tight_coupling = 0.d0
    y                = 0.d0
    dydx             = 0.d0

    open(27, file='vandb_pol.dat')
    open(28, file='phi_theta_pol.dat')
    open(29, file='pol_quant.dat')
    open(30, file='neut_quant.dat')


    ! Propagate each k-mode independently

    write(*,*) "Integrating Perturbations."

    do k = 1, n_k

       write(*,*) "k =", k
       
       k_current = ks(k)  ! Store k_current as a global module variable
       ck = c*k_current

       ! Initialize equation set for tight coupling
       y_tight_coupling(1)  = delta(1,k)
       y_tight_coupling(2)  = delta_b(1,k)
       y_tight_coupling(3)  = v(1,k)
       y_tight_coupling(4)  = v_b(1,k)
       y_tight_coupling(5)  = Phi(1,k)
       y_tight_coupling(6)  = Theta(1,0,k)
       y_tight_coupling(7)  = Theta(1,1,k)
       y_tight_coupling(13) = Theta_p(1,0,k)
       y_tight_coupling(14) = Theta_p(1,1,k)
       y_tight_coupling(20) = Nu(1,0,k)
       y_tight_coupling(21) = Nu(1,1,k)
       y_tight_coupling(22) = Nu(1,2,k)
       y_tight_coupling(23) = Nu(1,3,k)
       y_tight_coupling(24) = Nu(1,4,k)
       y_tight_coupling(25) = Nu(1,5,k)
       y_tight_coupling(26) = Nu(1,6,k)
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)
       
       ! write(*,*) x_tc
       
       ! Write initial values to file for k=1,10,30,50,80,100
       do i = 1,6
          if (k == prints(i)) then
             write(27,'(5(E17.8))') x_t(1), delta(1,k), delta_b(1,k), v(1,k), v_b(1,k)
             write(28,'(5(E17.8))') x_t(1), Phi(1,k), Psi(1,k), Theta(1,0,k), Theta(1,1,k)
             write(29,'(3(E17.8))') Theta_p(1,0,k), Theta_p(1,1,k), Theta_p(1,2,k)
             write(30,'(7(E17.8))') Nu(1,0,k), Nu(1,1,k), Nu(1,2,k), Nu(1,3,k), Nu(1,4,k), Nu(1,5,k), Nu(1,6,k)
          end if
       end do  


       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations       
       j=2
       do while (x_t(j) < x_tc)
            x     = x_t(j)
            a     = exp(x)
            bleta = get_eta(x)
            H     = get_H_p(x)
            ckH   = ck/H
            dtau  = get_dtau(x)

            write(*,*) "------ y_tc ------"

            do i=1,26
               write(*,*) i, y_tight_coupling(i)
            end do
             
            ! Solve next evolution step
            call odeint(y_tight_coupling,x_t(j-1),x,eps,h1,hmin,deriv_tc,bsstep,output)

             
            ! Save variables
            delta(j,k)     = y_tight_coupling(1)
            delta_b(j,k)   = y_tight_coupling(2)
            v(j,k)         = y_tight_coupling(3)
            v_b(j,k)       = y_tight_coupling(4)
            Phi(j,k)       = y_tight_coupling(5)
            Theta(j,0,k)   = y_tight_coupling(6)
            Theta(j,1,k)   = y_tight_coupling(7)
            Theta(j,2,k)   = -8.d0*ckH/(15.d0*dtau)*Theta(j,1,k)

            Theta_p(j,0,k) = y_tight_coupling(13)
            Theta_p(j,1,k) = y_tight_coupling(14)
            Theta_p(j,2,k) = 0.25d0*Theta(j,2,k)
            

            Psi(j,k)       = -Phi(j,k) - 12.d0*(H_0/(ck*a))**2.d0*Omega_r*Theta(j,2,k)

            Pi_c           = Theta(j,2,k) + Theta_p(j,0,k) + Theta_p(j,2,k)

            ! Task: Store derivatives that are required for C_l estimation

            call deriv_tc(x_t(j),y_tight_coupling,dydx)

            write(*,*) '------ dydx ------'
            do i =1,26
               write(*,*) i, dydx(i)
            end do

            dv_b(j,k)      = dydx(4)
            dPhi(j,k)      = dydx(5)
            
            dTheta(j,0,k)  = dydx(6)
            dTheta(j,1,k)  = dydx(7)
            dTheta(j,2,k)  = 2.d0/5.d0*ckH*Theta(j,1,k) -&
                             3.d0*ckH/(5.d0)*Theta(j,3,k)+dtau*(Theta(j,2,k)-0.1d0*Pi_c)

            dTheta_p(j,0,k)= dydx(13)
            dTheta_p(j,1,k)= dydx(14)
            dTheta_p(j,2,k)= 2.d0/5.d0*ckH*Theta_p(j,1,k) -&
                             3.d0*ckH/5.d0*Theta(j,3,k) + dtau*(Theta_p(j,2,k)-0.1d0*Pi_c)
           
            dNu(j,0,k)     = dydx(20)
            dNu(j,1,k)     = dydx(21)
            dNu(j,2,k)     = dydx(22)
            dNu(j,3,k)     = dydx(23)
            dNu(j,4,k)     = dydx(24)
            dNu(j,5,k)     = dydx(25)
            dNu(j,6,k)     = dydx(26)            

            dPsi(j,k)      = -dPhi(j,k) - 12.d0*H_0**2.d0/(ck*a)**2.d0*Omega_r*&
                             (-2.d0*Theta(j,2,k)+dTheta(j,2,k))

            ! Write values to file for k=1,10,30,50,80,100
            do i = 1,6
               if (k == prints(i)) then
                  write(27,'(5(E17.8))') x_t(j), delta(j,k), delta_b(j,k), v(j,k), v_b(j,k)
                  write(28,'(5(E17.8))') x_t(j), Phi(j,k), Psi(j,k), Theta(j,0,k), Theta(j,1,k)
                  write(29,'(3(E17.8))') Theta_p(j,0,k), Theta_p(j,1,k), Theta_p(j,2,k)
                  write(30,'(7(E17.8))') Nu(j,0,k),Nu(j,1,k),Nu(j,2,k),Nu(j,3,k),Nu(j,4,k),Nu(j,5,k),Nu(j,6,k)
               end if
            end do
            j = j+1
       end do  
       j_tc = j
        
       STOP

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       y(1:7) = y_tight_coupling(1:7)
       y(8)   = -20.d0*ckH/(45.d0*dtau)*Theta(i,1,k)

       
 
       do l = 3, lmax_int
          y(6+l) = -l*ckH/((2.d0*l+1.d0)*dtau)*y(6+l-1)
       end do

       do i = j_tc,n_t
          

             x = x_t(i)
             a = exp(x)
             bleta = get_eta(x)
             H   = get_H_p(x)
             ckH = ck/H
             dtau = get_dtau(x)
   
          ! Task: Integrate equations from tight coupling to today
          call odeint(y,x_t(i-1),x,eps,h1,hmin,deriv,bsstep,output)

          ! Task: Store variables at time step i in global variables
          delta(i,k)   = y(1) 
          delta_b(i,k) = y(2)
          v(i,k)       = y(3)
          v_b(i,k)     = y(4)
          Phi(i,k)     = y(5)

          do l = 0, lmax_int
             Theta(i,l,k) = y(6+l)
          end do

          Psi(i,k)     = - Phi(i,k) - 12.d0*(H_0/(ck*a))**2.d0*Omega_r*Theta(i,2,k)

          ! Task: Store derivatives that are required for C_l estimation
          call deriv(x_t(i),y,dydx)
          
          dPhi(i,k)     = dydx(5)
          dv_b(i,k)     = dydx(4)
          
          do l=0,lmax_int
             dTheta(i,:,k) = dydx(6+l) 
          end do

          dPsi(i,k)     = -dPhi(i,k) - (12.d0*H_0**2.d0)/(ck*a)**2.d0*&
                           Omega_r*(-2.d0*Theta(i,2,k)+dTheta(i,2,k))

          do j = 1,6
             if (k == prints(j)) then
                write(27,'(5(E17.8))') x_t(i), delta(i,k), delta_b(i,k), v(i,k), v_b(i,k)
                write(28,'(5(E17.8))') x_t(i), Phi(i,k), Psi(i,k), Theta(i,0,k), Theta(i,1,k)
             end if
          end do

       end do
    end do

    close(27)
    close(28)

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns

  ! -------------------- derivative subroutines ------------------------

  subroutine deriv_tc(x,y_tc,dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y_tc
    real(dp), dimension(:), intent(out) :: dydx
    
    real(dp) :: d_delta
    real(dp) :: d_delta_b
    real(dp) :: d_v
    real(dp) :: q,R

    real(dp) :: delta,delta_b,v,v_b,Phi,Theta0,Theta1,Theta2
    real(dp) :: Psi,dPhi,dTheta0,dv_b,dTheta1
    real(dp) :: dThetap0, dThetap1,dN0,dN1,dN2,dN3,dN4,dN5,dN6
    real(dp) :: dtau, ddtau,a,H_p,dH_p,ckH_p,bleta
    real(dp) :: Thetap0,Thetap1,Thetap2,N0,N1,N2,N3,N4,N5,N6,Pi_d
    
    delta     = y_tc(1)
    delta_b   = y_tc(2)
    v         = y_tc(3)
    v_b       = y_tc(4)
    Phi       = y_tc(5)
    Theta0    = y_tc(6)
    Theta1    = y_tc(7)
    Thetap0   = y_tc(13)
    Thetap1   = y_tc(14)
    N0        = y_tc(20)
    N1        = y_tc(21)
    N2        = y_tc(22)
    N3        = y_tc(23)
    N4        = y_tc(24)
    N5        = y_tc(25)
    N6        = y_tc(26)

    dtau      = get_dtau(x)
    ddtau     = get_ddtau(x)
    a         = exp(x)
    H_p       = get_H_p(x)
    dH_p      = get_dH_p(x)
    ckH_p     = c*k_current/H_p
    bleta     = get_eta(x)
    
    Theta2    = -8.d0*ckH_p/(15.d0*dtau)*Theta1
    Thetap2   = 0.25d0*Theta2

    R         = (4.d0*Omega_r)/(3.d0*Omega_b*a)

    Pi_d      = Thetap0+ Theta2 + Thetap2

    Psi       = -Phi - 12.d0*(H_0/(c*k_current*a))**2.d0*(Omega_r*Theta2+Omega_nu*N2)

    dPhi      = Psi - (ckH_p**2.d0)/3.d0*Phi + (H_0/H_p)**2.d0/2.d0*(Omega_m/a*delta + &
               Omega_b/a*delta_b + 4.d0*Omega_r*Theta0/a**2.d0 + 4.d0*Omega_nu/a**2.d0*N0)

    dTheta0   = -ckH_p*Theta1 - dPhi
    
    d_delta   = ckH_p*v - 3.d0*dPhi  

    d_delta_b = ckH_p*v_b - 3.d0*dPhi

    d_v       = -v - ckH_p*Psi

    q         = (-((1.d0-2.d0*R)*dtau + (1.d0+R)*ddtau)*(3.d0*Theta1 + v_b) - ckH_p*Psi + &
                (1.d0-dH_p/H_p)*ckH_p*(-Theta0-2.d0*Theta2)&
                - ckH_p*dTheta0)/((1.d0+R)*dtau+dH_p/H_p - 1.d0)

    dv_b      = (1.d0/(1.d0+R))*(-v_b-ckH_p*Psi + R*(q+ckH_p*(-Theta0+2.d0*Theta2)-ckH_p*Psi))

    dTheta1   = (1.d0/3.d0)*(q-dv_b)
    
    dThetap0  = -ckH_p*Thetap1 + dtau*(Thetap0-0.5d0*Pi_d)
    dThetap1  = ckH_p/3.d0*Thetap0 + dtau*Thetap1 - 2.d0*ckH_p/3.d0*Thetap2

    dN0       = -ckH_p*N1 - dPhi
    dN1       = ckH_p*N0 - 2.d0*ckH_p/3.d0*N2 + ckH_p*Psi/3.d0
    dN2       = (2.d0*ckH_p)/(5.d0)*N1 - (3.d0)*ckH_p/(5.d0)*N3
    dN3       = (3.d0*ckH_p)/(7.d0)*N2 - (4.d0)*ckH_p/(7.d0)*N4
    dN4       = (4.d0*ckH_p)/(9.d0)*N3 - (5.d0)*ckH_p/(9.d0)*N5
    dN5       = (5.d0*ckH_p)/(11.d0)*N4 - (6.d0)*ckH_p/(11.d0)*N6
    dN6       = (ckH_p)/(7.d0)*N2 - c*N6*7.d0/(H_p*bleta)    

    ! Output
    dydx(1)   = d_delta
    dydx(2)   = d_delta_b
    dydx(3)   = d_v
    dydx(4)   = dv_b
    dydx(5)   = dPhi
    dydx(6)   = dTheta0
    dydx(7)   = dTheta1
    dydx(13)  = dThetap0
    dydx(14)  = dThetap1
    dydx(20)  = dN0
    dydx(21)  = dN1
    dydx(22)  = dN2
    dydx(23)  = dN3
    dydx(24)  = dN4
    dydx(25)  = dN5
    dydx(26)  = dN6

  end subroutine deriv_tc


  subroutine deriv(x,y,dydx)
    use healpix_types
    implicit none
    
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    real(dp) :: d_delta
    real(dp) :: d_delta_b
    real(dp) :: d_v
    real(dp) :: R

    real(dp) :: delta,delta_b,v,v_b,Phi,Theta0,Theta1,Theta2,Theta3,Theta4,Theta5,Theta6
    real(dp) :: Psi, dPhi,dTheta0,dv_b,dTheta1,dTheta2
    real(dp) :: a,H_p,ckH_p,dtau,bleta
    integer(i4b) :: l

    delta     = y(1)
    delta_b   = y(2)
    v         = y(3)
    v_b       = y(4)
    Phi       = y(5)
    Theta0    = y(6)
    Theta1    = y(7)
    Theta2    = y(8)
    Theta3    = y(9)
    Theta4    = y(10)
    Theta5    = y(11)
    Theta6    = y(12)

    a         = exp(x)
    H_p       = get_H_p(x)
    ckH_p     = c*k_current/H_p
    dtau      = get_dtau(x)
    bleta     = get_eta(x)

    R         = (4.d0*Omega_r)/(3.d0*Omega_b*a)
    Psi       = -Phi - 12.d0*H_0*H_0/((c*k_current*a)*(c*k_current*a))*Omega_r*Theta2

    dPhi      = Psi - (ckH_p**2.d0)/3.d0*Phi + H_0**2.d0/(2.d0*H_p**2.d0)*(Omega_m/a*delta + &
                Omega_b/a*delta_b + 4.d0*Omega_r*Theta0/a**2.d0)

    dTheta0   = -ckH_p*Theta1 - dPhi
    d_delta   = ckH_p*v   - 3.d0*dPhi  
    d_delta_b = ckH_p*v_b - 3.d0*dPhi
    d_v       = -v - ckH_p*Psi
    
    dv_b      = -v_b - ckH_p*Psi + dtau*R*(3.d0*Theta1 + v_b)

    dTheta1   = ckH_p*Theta0/3.d0 - 2*ckH_p*Theta2/3.d0 + ckH_p*Psi/3.d0 + dtau*(Theta1 + v_b/3.d0)
    dTheta2   = 2.d0/5.d0*ckH_p*Theta1 - 3.d0/5.d0*ckH_p*Theta3+dtau*0.9d0*Theta2

    do l=3,lmax_int-1
       dydx(6+l) = l/(2.d0*l+1.d0)*ckH_p*y(5+l) - (l+1.d0)/(2.d0*l+1.d0)*ckH_p*y(7+l) + dtau*y(6+l)
    end do 

    dydx(6+lmax_int) = ckH_p*y(5+lmax_int) - c*(lmax_int+1.d0)/(H_p*bleta)*y(6+lmax_int) &
                       + dtau*y(6+lmax_int)

    ! Output
    dydx(1)   = d_delta
    dydx(2)   = d_delta_b
    dydx(3)   = d_v
    dydx(4)   = dv_b
    dydx(5)   = dPhi
    dydx(6)   = dTheta0
    dydx(7)   = dTheta1
    dydx(8)   = dTheta2

  end subroutine deriv

  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x < x(start of recombination)

  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    integer(i4b)          :: i,n
    real(dp)              :: x
    n = 1d4
    do i=0,n
       x = x_init + i*(0.d0-x_init)/n
       if (x < x_start_rec .and. abs(c*k/(get_H_p(x)*get_dtau(x))) < 0.1d0 .and. &
            abs(get_dtau(x)) > 10.d0) then
          get_tight_coupling_time = x
       end if
    end do
  end function get_tight_coupling_time

end module evolution_mod_w_pol
