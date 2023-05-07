!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     int.f90: 
!     integrates f(x)=exp(-x**2) in the interval [0,1]
!     using Monte Carlo importance sampling algorithm
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module intmod
  public :: f, f_, mc
contains

  ! function to be integrated
  !
  function f(x)
    implicit none
    real :: f
    real, intent(in) :: x
    f = exp(-x**2)
    return
  end function f
  
  FUNCTION rand_exp(i)

	REAL :: rand_exp
	REAL :: r
	REAL, intent(in)  :: i
	26 continue
	DO
		CALL RANDOM_NUMBER(r)
		IF (r > 0) EXIT
	END DO

	rand_exp = -((1/i)*LOG(r))
	IF (rand_exp > 1.) go to 26
	
	RETURN
  END FUNCTION rand_exp

  
  ! function to be integrated
  !
  function f_(x)
    implicit none
    real :: f_
    real, intent(in) :: x
    f_ = (exp(-x**2.)/exp(-x))*(1.-exp(-1.))
    return
  end function f_

  ! Monte Carlo
  !
  function mc(i, distribution)
    implicit none
	real, dimension(2) :: mc
    integer, intent(in) :: i
	character (len = *),intent(in) :: distribution
    integer :: n
	real :: r
    mc(1) = 0.
	mc(2) = 0.
	if (distribution=='Uniform') then
		do n = 1,i
			call random_number(r)
			mc(1)=f(r)/i+mc(1)
			mc(2)=f(r)**2/i+mc(2)
		end do
	else if(distribution=='Exponential') then
		do n = 1,i
		    r=rand_exp(1.)
			mc(1)=f_(r)/i+mc(1)
			mc(2)=f_(r)**2/i+mc(2)
		end do
	end if
	return
  end function mc

end module intmod

program int
  use intmod
  implicit none
  integer :: sizer
  integer, dimension(:), allocatable :: seed
  real, dimension(2) :: r , r_
  real :: theo=0., t0, t1,PI
  real :: var, var_, std, std_, error, error_
  integer :: i, n
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  seed=[1,2,3,4,5,6,7,8]
  call random_seed(put=seed)

 
  PI=4.D0*DATAN(1.D0)
  theo = erf(1.)*sqrt(PI)*0.5
  
  print*,'Exact Value =',theo
  
  open(unit=7,file='Uniform.dat',status='unknown')
  open(unit=8,file='Exponential.dat',status='unknown')
  
  write(7,*)"N, Integral(Uniform), Variance(Uniform), STD(Uniform), Error(Uniform)"
  write(8,*)"N, Integral(Exponential), Variance(Exponential), STD(Exponential), Error(Exponential)"
  
  call cpu_time(t0)
  do i = 10,10000
     r = mc(i,'Uniform')
	 var=r(2)-r(1)**2
	 std=sqrt(var/i)
     error = abs(r(1)-theo)
	 
	 r_ = mc(i,'Exponential')
	 var_=r_(2)-r_(1)**2
	 std_=sqrt(var_/i)
     error_ = abs(r_(1)-theo)
	 
     write(7,'(i5,4(5x,f10.6),4(7x,f10.6),4(5x,f10.6),4(5x,f10.6))') i, r(1), var, std, error
	 write(8,'(i5,4(5x,f10.6),4(7x,f10.6),4(5x,f10.6),4(5x,f10.6))') i, r_(1), var_, std_, error_
  end do
  call cpu_time(t1)
  print*,"Calculated with Uniform sampling:", r(1)
  print*,"Calculated with Exponential sampling:", r_(1)
  print*,"Total Time spent:",t1-t0
  close(7)
  close(8)

  print*,'Data saved.'
  stop

end program int
