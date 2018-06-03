module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none


  real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0


  integer(i4b)                           :: n_t, n1, n2,n3     ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values
  real(dp),    allocatable, dimension(:) :: eta_t              ! Grid of relevenat eta-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: z_eta              ! Grid points for z_eta
  real(dp),    allocatable, dimension(:) :: a_eta              ! Grid points for a_eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

  real(dp)                               :: rho_crit0          ! critical density today
  real(dp)                               :: rho_m0             ! matter density today
  real(dp)                               :: rho_b0             ! baryon density today
  real(dp)                               :: rho_r0             ! radiation density today
  real(dp)                               :: rho_nu0            ! neutrino density today
  real(dp)                               :: rho_lambda0        ! vacuum energy density today

  real(dp),    allocatable, dimension(:) :: rho_m              ! matter density grid
  real(dp),    allocatable, dimension(:) :: rho_b              ! baryon density grid
  real(dp),    allocatable, dimension(:) :: rho_r              ! radiation density grid
  real(dp),    allocatable, dimension(:) :: rho_nu             ! neutrino density grid
  real(dp),    allocatable, dimension(:) :: rho_lambda         ! vacuum energy density grid

  real(dp),    allocatable, dimension(:) :: Omega_mx           ! Relative densities
  real(dp),    allocatable, dimension(:) :: Omega_bx           ! ''
  real(dp),    allocatable, dimension(:) :: Omega_rx           ! ''
  real(dp),    allocatable, dimension(:) :: Omega_nux          ! ''
  real(dp),    allocatable, dimension(:) :: Omega_lambdax      ! ''

  real(dp),    allocatable, dimension(:) :: Hx                 ! Hubble constant as a func of x

  real(dp)                               :: m2mpc              ! Value for changing m to Mpc
  real(dp)                               :: kmmpc


contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n
    real(dp)     :: dx,x_eta1, x_eta2, a_init,a_end,eta_init
    real(dp)     :: eps,hmin,yp1,ypn,h1
    real(dp),    allocatable, dimension(:) :: dydx

    ! Define two epochs, 1) during and 2) after recombination.

    !-----------------------------------------------------------------------------------
    ! Initialize the values for variables that are used during the calculations

    n1          = 1000                      ! Number of grid points before recombination
    n2          = 200                       ! Number of grid points during recombination
    n3          = 300                       ! Number of grid points after recombination
    n_t         = n1+n2+n3                  ! Total number of grid points

    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today

    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-8                     ! Start value of a for eta evaluation
    a_end       = 1.d0
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    eps  = 1.d-10
    hmin = 0.d0

    m2mpc = 3.2408d-23                      ! 3.086*10^22m = 1Mpc
    kmmpc = 3.086d19                        ! ''

    yp1   = 1.d30
    ypn   = 1.d30


    !-----------------------------------------------------------------------------------
    ! Task: Fill in x and a grids
    allocate(x_t(n_t))
    allocate(a_t(n_t))
    allocate(eta(n_eta))
    allocate(x_eta(n_eta))
    allocate(a_eta(n_eta))
    allocate(z_eta(n_eta))
    allocate(eta2(n_eta))

    allocate(rho_m(n_eta))
    allocate(rho_b(n_eta))
    allocate(rho_r(n_eta))
    allocate(rho_nu(n_eta))
    allocate(rho_lambda(n_eta))

    allocate(Omega_mx(n_eta))
    allocate(Omega_bx(n_eta))
    allocate(Omega_rx(n_eta))
    allocate(Omega_nux(n_eta))
    allocate(Omega_lambdax(n_eta))
    allocate(Hx(n_eta))


    do i = 1,n1                           ! Fill recombination interval
       x_t(i) = x_eta1 + (i-1)*(x_start_rec-log(a_init))/(n1-1)
    end do

    do i = 1,n2                             ! Fill post recombination interval
       x_t(n1+i) = x_start_rec + (i)*(x_end_rec-x_start_rec)/(n2)
    end do

    do i = 1,n3
       x_t(n1+n2+i) = x_end_rec + i*(-x_end_rec)/(n3)
    end do

    a_t = exp(x_t)                          ! Fill the a_t-grid using the x_t values


    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    
    ! Compute the x_eta values
    x_eta(1) = x_eta1
    do i = 1,n_eta-1
       x_eta(i+1) = x_eta1 + i*(x_eta2-x_eta1)/(n_eta-1)
    end do

    eta(1) = c*1.d-10/(H_0*sqrt(Omega_r+Omega_nu)) ! early universe is rad/neutrino dominated

    h1 = abs(1.d-2*(x_eta(1)-x_eta(2)))

    allocate(dydx(1))

    do i =2,n_eta
       eta(i) = eta(i-1)
       call odeint(eta(i:i),x_eta(i-1),x_eta(i),eps,h1,hmin,eta_derivs,bsstep,output)
    end do
    
    call spline(x_eta,eta,yp1,ypn,eta2)


    ! Compute the a_eta and z_eta values
    a_eta = exp(x_eta)
    z_eta = 1.d0/a_eta - 1.d0

    !--------------------------------------------------------------------------------------
    ! Density calculations 
    
    rho_crit0  = 3.d0*H_0**2.d0/(8.d0*pi*G_grav)
    rho_m0     = Omega_m      * rho_crit0
    rho_b0     = Omega_b      * rho_crit0
    rho_r0     = Omega_r      * rho_crit0
    rho_nu0    = Omega_nu     * rho_crit0
    rho_lambda0= Omega_lambda * rho_crit0

    do i=1,n_eta
       Hx(i) = get_H(x_eta(i))*kmmpc
       Omega_mx(i)      = Omega_m     *H_0**2.d0/Hx(i)**2.d0 * a_eta(i)**-3.d0
       Omega_bx(i)      = Omega_b     *H_0**2.d0/Hx(i)**2.d0 * a_eta(i)**-3.d0
       Omega_rx(i)      = Omega_r     *H_0**2.d0/Hx(i)**2.d0 * a_eta(i)**-4.d0
       Omega_nux(i)     = Omega_nu    *H_0**2.d0/Hx(i)**2.d0 * a_eta(i)**-4.d0
       Omega_lambdax(i) = Omega_lambda*H_0**2.d0/Hx(i)**2.d0
    end do

    ! Print the density values to file

    open(50, file='densities.dat')
      do i=1,n_eta
         write(50,'(5(E17.8))') Hx(i),Omega_mx(i),Omega_bx(i),Omega_rx(i),Omega_lambdax(i)
      end do
    close(50)



    !-------------------------------------------------------------------------------------
    ! Calculate the eta values, initializing eta(1) using the assumption that the universe
    ! is radiation dominated before recomination (rad/matter equality).


    ! write(*,*) eta ! eta check

    ! Print x,a, and eta values to file for plotting
    open(54, file='data.dat')
       do i=1,n_eta
          write(54,'(4(E17.8))') x_eta(i),a_eta(i),eta(i),z_eta(i)
       end do
    close(54)



  end subroutine initialize_time_mod

  !-----------------------------------------------------------------------------------

  subroutine eta_derivs(x,eta,dydx) ! Defining the derivative d(eta)/da

  use healpix_types
     implicit none

     real(dp),               intent(in)  :: x
     real(dp), dimension(:), intent(in)  :: eta
     real(dp), dimension(:), intent(out) :: dydx
     real(dp)                            :: H_p

     H_p = get_H_p(x)
     dydx = c/(H_p)

  end subroutine eta_derivs


  subroutine output(x,y)
     use healpix_types
     implicit none

     real(dp),               intent(in)  :: x
     real(dp), dimension(:), intent(in)  :: y

  end subroutine output 

!-------------------------------Done with subroutines----------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!---------------------------------Define Functions-------------------------------------------------

  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
    real(dp)             :: a
    a = exp(x)

    get_H = H_0*sqrt((Omega_b+Omega_m)*a**-3.d0 + (Omega_r + Omega_nu)*a**-4.d0 + Omega_lambda)

  end function get_H


  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    real(dp)             :: a
    a = exp(x)

    get_H_p = a*get_H(x)

  end function get_H_p


  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
    get_dH_p = H_0/2.d0*1/sqrt((Omega_m+Omega_b)*exp(-x)+Omega_r*exp(-2.d0*x) &
    + Omega_lambda*exp(2.d0*x)) * (-(Omega_m+Omega_b)*exp(-x)-2.d0*Omega_r*exp(-2.d0*x) &
    + 2.d0*Omega_lambda*exp(2.d0*x))

  end function get_dH_p


  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none
    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    get_eta = splint(x_eta,eta,eta2,x_in)
  end function get_eta

end module time_mod
