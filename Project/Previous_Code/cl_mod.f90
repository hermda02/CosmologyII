module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

  real(dp),     pointer,     dimension(:)       :: x_hires, k_hires, l_hires, cl_hires
  real(dp),     pointer,     dimension(:,:)     :: S, S2, j_l, j_l2
  real(dp),     allocatable, dimension(:)       :: bessel,z_spline
  real(dp),     allocatable, dimension(:,:)     :: Theta_l,integ2,integ
  integer(i4b), allocatable, dimension(:)       :: ls
  
contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    real(dp)     :: k_min, k_max, dk, l_max
    real(dp),     allocatable, dimension(:)       :: integrand,int_test,int2_test
    real(dp),     pointer,     dimension(:)       :: x_arg, cls, cls2, ls_dp
    real(dp),     pointer,     dimension(:)       :: k, x, cl_int
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     allocatable, dimension(:)       :: j_l_spline, j_l_spline2

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    n_spline = 5400
    l_num = 44

    ! Allocate all necessary arrays for C_l estimation
    allocate(S(n_x_hires,n_k_hires))
    allocate(x_hires(n_x_hires),k_hires(n_k_hires))

    allocate(ls(l_num))
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, dummy argument of j_l(z)

    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    allocate(Theta_l(l_num,n_k_hires))
    allocate(integ(l_num,n_x_hires))
    allocate(integ2(l_num,n_k_hires))

    allocate(cl_int(l_num))
    allocate(cls(l_num))
    allocate(cls2(l_num))

    allocate(ls_dp(l_num))

    allocate(int_test(n_x_hires))
    allocate(int2_test(n_k_hires))


    ! Open C_l files to write to:
    !----------------------------
    open(11, file='integrand_l_100.dat')
    open(10, file='thetal_l_100.dat')
    open(12, file='hi_res_C_ls.dat')

    ! Set up which l's to compute
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Convert the l's to double precision, and create high-resolution l arrays
    do l=1,l_num
       ls_dp(l) = ls(l)
    end do

    l_max = int(maxval(ls))

    allocate(l_hires(int(l_max)))
    allocate(cl_hires(int(l_max)))

    do l=1,l_max
       l_hires(l) = l
    end do

    ! Task: Get source function from evolution_mod

    call get_hires_source_function(x_hires,k_hires,S)

    ! Initialize values for trapezoidal integration
    x_min = x_hires(1)
    x_max = x_hires(n_x_hires)
    dx    = (x_max-x_min)/n_x_hires    
    k_min = k_hires(1)
    k_max = k_hires(n_k_hires)
    dk    = (k_max-k_min)/n_k_hires    
 
    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.

    do i=1,n_spline
       z_spline(i) = 0.d0 + (i-1)*3500.d0/(n_spline-1.d0)
    end do

    inquire(file='bessel.unf', exist=exist)
    if (exist) then
       write(*,*) 'Bessel file found.'
       open(13,file='bessel.unf', form='unformatted')
           read(13) j_l
       close(13)
    else
       write(*,*) 'Calculating Bessel function values.'
       do i=2,n_spline
          do l=1,l_num
             j_l(1,l) = 0.d0
             call sphbes(ls(l),z_spline(i),j_l(i,l))
          end do
       end do
    
       open(13,file='bessel.unf', form='unformatted')
           write(13) j_l
       close(13)
    end if

    do l=1,l_num
       call spline(z_spline,j_l(:,l),yp1,ypn,j_l2(:,l))
    end do

    !write(*,*) splint(z_spline,j_l(:,2), j_l2(:,2), 0.2*H_0/c*get_eta(0.d0))
    !write(*,*) get_eta(0.d0), 0.2*H_0/c*get_eta(0.d0)
   
    !| Overall task: Compute the C_l's for each given l
    !|--------------------------------------------------

    !|--------------------------------------
    !| Integrand test for l=100, k = 159.988
    !|--------------------------------------
    do i=1,n_x_hires
       int_test(i) = S(i,2000)*j_lfunc(17,k_hires(2000),x_hires(i))
    end do

    write(*,*) 'Computing C_ls!'

    do l = 1,l_num
       write(*,*) 'l = ',ls(l)
       
       ! Task: Compute the transfer function, Theta_l(k)
       ! We will compute the integral by using the trapezoidal method,
       ! summing over all x-values of the Bessel weighted Source Function

       do j=1,n_k_hires

          do i=1,n_x_hires
             integ(l,i) = S(i,j)*j_lfunc(l,k_hires(j),x_hires(i))
             Theta_l(l,j) = Theta_l(l,j) + integ(l,i)
          end do
          Theta_l(l,j) = dx*Theta_l(l,j)

          ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
          !--------------------------------------------------------------------------
          integ2(l,j) = (c*k_hires(j)/H_0)**(n_s-1)*Theta_l(l,j)**(2.d0)/k_hires(j)
          cl_int(l) = cl_int(l) + integ2(l,j)
       end do

       ! Task: Store C_l in an array. Optionally output to file
       !-------------------------------------------------------
       cls(l) = dk*cl_int(l)*ls(l)*(ls(l)+1.d0)/(2.d0*pi)

       write(*,*) 'C_l = ', cls(l)
    end do

    do i=1,n_k_hires
       int2_test(i) = (c*k_hires(i)/H_0)**(n_s-1)*Theta_l(17,i)**(2.d0)/k_hires(i)
    end do


    !Integrand testing apparatus:
!----------------------------------------------------------------
                                                             !|
    do j=1,n_k_hires                                         !|
       if (abs(int2_test(j)) < 1.d-99) then                  !|
          int2_test(j) = 0.d0                                !|
       endif                                                 !|
       if (abs(int_test(j)) < 1.d-99) then                   !|
          int_test(j) = 0.d0                                 !|
       endif                                                 !|
       if (abs(S(j,2000)) < 1.d-99) then                     !|
          S(j,2000) = 0.d0                                   !|
       endif                                                 !|
       write(10,'(2(E18.7))') c*k_hires(j)/H_0, int2_test(j) !|
       write(11,'(2(E18.7))') x_hires(j), int_test(j)        !|
    end do                                                   !|
                                                             !|
    close(10)                                                !|
    close(11)                                                !|
!----------------------------------------------------------------


    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l

    ! Now spline the C_l's
    write(*,*) 'Splining C_ls'
    call spline(ls_dp,cls,yp1,ypn,cls2)
    write(*,*) 'C_l spline complete!'

    ! Create high-resolution C_l array
    do l=2,l_max
       cl_hires(l) = splint(ls_dp,cls,cls2,l_hires(l))
       ! write(*,*) 'l = ', l_hires(l), 'C_l = ', cl_hires(l)
       write(12,'(2(E17.8))') l_hires(l), cl_hires(l)
    end do

    close(12)    

  end subroutine compute_cls
  

  function j_lfunc(l,k,x)
    implicit none
    integer(i4b),  intent(in) :: l
    real(dp),      intent(in) :: x,k
    real(dp)                  :: j_lfunc

    j_lfunc = splint(z_spline,j_l(:,l),j_l2(:,l),k*(get_eta(0.d0)-get_eta(x)))

  end function j_lfunc


end module cl_mod
