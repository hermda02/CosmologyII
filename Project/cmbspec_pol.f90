program cmbspec_pol
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod_w_pol
  !use cl_mod
  implicit none

  ! Initialize time grids
  call initialize_time_mod
  !write(*,*) "Hello World! :)"
  ! Output to file desired quantities here
  call initialize_rec_mod

  call initialize_perturbation_eqns
  call integrate_perturbation_eqns

  ! call compute_cls

end program cmbspec_pol