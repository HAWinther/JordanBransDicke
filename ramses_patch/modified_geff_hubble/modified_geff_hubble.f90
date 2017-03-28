module modified_geff_hubble
  implicit none
  public :: GeffOverG, HubbleFunction, init_modified_geff_hubble_simulation

contains

  ! This function is called everywhere we have a 'fourpi = 1.5 Omega_m a' factor
  ! This denotes the prefactor of the Poisson equation so we include GeffG mod
  ! by taking fourpi => fourpi * GeffOverG(a) 
  function GeffOverG(a)
    use class_QubicSpline
    use SplineContainer
    implicit none
    real(kind=8) :: GeffOverG, a

#if defined(MODIFIEDGEFFHUBBLE)
    ! The spline contains GeffG( log(a) )
    GeffOverG = LookupSpline(GeffSpline, log(a))
#else
    GeffOverG = 1.0d0
#endif
    return
  end function
  
  ! This function is only called from init_time.f90
  ! and not called at all if MODIFIEDGEFFHUBBLE is not defined
  function HubbleFunction(a)
    use amr_parameters, only: omega_m, omega_k, omega_l
    use class_QubicSpline
    use SplineContainer
    implicit none
    real(kind=8) :: HubbleFunction, a
    
#if defined(MODIFIEDGEFFHUBBLE)
    ! The spline contains log[ H( log (a) ) ]
    HubbleFunction = exp(LookupSpline(HubbleSpline, log(a)))
#else
    HubbleFunction = sqrt(omega_m/a**3 + omega_k/a**2 + omega_l)
#endif
    return
  end function
 
  ! Reads the file in namelist/amr/filename_geff_hubble_data
  ! whose format is as follows:
  !-------------------------------------
  ! [One line header]
  ! n    
  ! log(a_1)   GeffG(a_1)   H(a_1)/H0
  ! log(a_2)   GeffG(a_2)   H(a_2)/H0
  ! ...
  ! log(a_n)   GeffG(a_n)   H(a_n)/H0
  !-------------------------------------
  ! where we require a_(i+1) > a_i
  ! Then splines the two functions and makes the availiable from
  ! calling the two functions above.
  subroutine init_modified_geff_hubble_simulation
    use class_QubicSpline
    use SplineContainer
    use amr_parameters, only: filename_geff_hubble_data, omega_m, omega_l, omega_k, aexp
    use amr_commons, only: myid
    implicit none
    integer :: ilun = 10, i, n, stat
    real(kind=8), dimension(:), allocatable :: x, Geffofx, logHofx
    real(kind=8) :: GeffGin, Hofxin
    real(kind=8) :: dx, zero
    character(len=1000) header
    
    ! Open the file containg GeffG(a) and H(a)/H0 
    open (unit=ilun, file=filename_geff_hubble_data, status='old', access='sequential', form='formatted', action='read', iostat=stat)
    if(stat .ne. 0)then
      write(*,*) 'Error: file not found ', filename_geff_hubble_data
      call clean_stop
      stop
    endif
    
    ! Read header from file and get number of points
    read(ilun, '(A)') header
    read(ilun, *) n

    ! Allocate temp-arrays
    allocate(x(1:n))
    allocate(Geffofx(1:n))
    allocate(logHofx(1:n))
   
    ! Read data from file
    do i = 1, n
      read(ilun, *) x(i), GeffGin, Hofxin
      Geffofx(i) = GeffGin
      logHofx(i) = log(Hofxin)
    enddo
    close(ilun)
    
    ! Sanity checks
    if(exp(x(1)) > aexp) then
      write(*,*) 'Error: the file with H(a), GeffG(a) has amin > aini : ', exp(x(1)), aexp
      call clean_stop
      stop
    endif
    do i = 1, n-1
      if(x(i) >= x(i+1) ) then
        write(*,*) 'Error: the log(a)-values in the file with H(a), GeffG(a) is not stricktly increasing'
        write(*,*) 'i      = ', i
        write(*,*) 'x(i)   = ', x(i)
        write(*,*) 'x(i+1) = ', x(i+1)
      endif
    enddo

    ! Create splines
    call CreateSpline(HubbleSpline, x, logHofx, n)
    call CreateSpline(GeffSpline, x, Geffofx, n)
  
    if(myid==1)then
      write(*,*) ''
      write(*,*) '==================================================================='
      write(*,*) 'Running a simulations with a modified Geff/G and/or Hubble equation' 
      write(*,*) '==================================================================='
      write(*,*) 'Reading data from : ', trim(filename_geff_hubble_data)
      write(*,*) 'Header            : ', trim(header)
      write(*,*) 'n_in_file         : ', n
      write(*,*) 'GeffG(a=0.01)     : ', GeffOverG(0.01d0)
      write(*,*) 'GeffG(a=0.1)      : ', GeffOverG(0.10d0)
      write(*,*) 'GeffG(a=1)        : ', GeffOverG(1.00d0)
      write(*,*) 'H(a=0.01) / Hlcdm : ', HubbleFunction(0.01d0) / sqrt(omega_m/0.01**3 + omega_k/0.01**2 + omega_l)
      write(*,*) 'H(a=0.1)  / Hlcdm : ', HubbleFunction(0.10d0) / sqrt(omega_m/0.10**3 + omega_k/0.10**2 + omega_l)
      write(*,*) 'H(a=1.0)  / Hlcdm : ', HubbleFunction(1.00d0) / sqrt(omega_m         + omega_k         + omega_l)
      if(HubbleSpline%linear_x_array) write(*,*) 'Linear x-array. Using direct lookup in splines'
      write(*,*) ''
    endif
  
    ! Free up memory
    deallocate(x)
    deallocate(Geffofx)
    deallocate(logHofx)
  end subroutine init_modified_geff_hubble_simulation
  
end module modified_geff_hubble
