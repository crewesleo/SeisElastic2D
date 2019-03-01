
program test_adjoint_sources
 
use constants
implicit none

 ! synthetics
 integer, parameter :: NSTEP = 1000
 real(kind=CUSTOM_REAL), parameter :: deltat=0.1
 real(kind=CUSTOM_REAL), parameter :: sigma=3.5
 real(kind=CUSTOM_REAL), parameter :: s0=50, d0=52
 real(kind=CUSTOM_REAL), parameter :: fc=100
 real(kind=CUSTOM_REAL), dimension(NSTEP) :: s,d,se,de,time

 ! window
 integer :: ntstart=301
 integer :: ntend=701
 real(kind=CUSTOM_REAL), parameter :: t0=0
 real(kind=CUSTOM_REAL), parameter :: f0=fc
 integer, parameter :: window_type=3   ! 1 -- boxcar; 2 -- welch ; 3 -- cosine


 ! misfit
 character(len=2) :: measurement_type(5)=['CC','WD','IP','ED','MT']
 real(kind=CUSTOM_REAL) :: misfit_output

 ! adjoint 
 logical :: compute_adjoint = .true.
 real(kind=CUSTOM_REAL), dimension(NSTEP) :: adj
 
 ! index
 integer :: i,itype


 !! gaussian functions
  ! time 
  do i=1,NSTEP
     time(i) = (i-1)*deltat
  enddo

  ! s
  call gaussmf(time,sigma,s0,NSTEP,s)
!  call gauspuls(time,NSTEP,fc,sigma,c,se,s )

  ! d 
  call gaussmf(time,sigma,d0,NSTEP,d)
!  call gauspuls(time,NSTEP,fc,sigma,c,de,d )


do itype=1,size(measurement_type)

!! initialization 
adj(:) = 0.d0 

!! misfit and adjoint source
call misfit_adj_AD(measurement_type(itype),d,s,NSTEP,deltat,f0,ntstart,ntend,&
                   window_type,compute_adjoint, &
                   misfit_output,adj)

print*, 'measurement_', trim(measurement_type(itype)) ,' = ', misfit_output
print*, 'squared misfit_', trim(measurement_type(itype)) ,' = ', misfit_output**2
print*

open(IOUT,file=trim(output_dir)//'/adj_'//trim(measurement_type(itype)),status='unknown')
do  i =  ntstart,ntend
    write(IOUT,*) time(i),d(i),s(i),adj(i)
enddo
close(IOUT)
enddo

end program test_adjoint_sources
