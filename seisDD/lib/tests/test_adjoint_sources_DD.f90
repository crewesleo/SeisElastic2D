
program test_adjoint_sources_DD
 
use constants
implicit none

 ! synthetics
 integer, parameter :: NSTEP = 1000
 real(kind=CUSTOM_REAL), parameter :: deltat=0.1
 real(kind=CUSTOM_REAL), parameter :: sigma=3.5
 real(kind=CUSTOM_REAL), parameter :: s10=50, d10=52, s20=51, d20=52
 real(kind=CUSTOM_REAL), parameter :: fc=100
 real(kind=CUSTOM_REAL), dimension(NSTEP) :: s1,d1,s2,d2,se,de,time

 ! window
 integer :: ntstart1=301
 integer :: ntend1=701
 integer :: ntstart2=301
 integer :: ntend2=701

 real(kind=CUSTOM_REAL), parameter :: t0=0
 real(kind=CUSTOM_REAL), parameter :: f0=fc
 integer, parameter :: window_type=3   ! 1 -- boxcar; 2 -- welch ; 3 -- cosine


 ! misfit
 character(len=2) :: measurement_type(1)=['CC']
 real(kind=CUSTOM_REAL) :: misfit_output

 ! adjoint 
 logical :: compute_adjoint = .true.
 real(kind=CUSTOM_REAL), dimension(NSTEP) :: adj1,adj2
 
 ! index
 integer :: i,itype


 !! gaussian functions
  ! time 
  do i=1,NSTEP
     time(i) = (i-1)*deltat
  enddo

  ! s1
  call gaussmf(time,sigma,s10,NSTEP,s1)

  ! d1
  call gaussmf(time,sigma,d10,NSTEP,d1)

  ! s2
  call gaussmf(time,sigma,s20,NSTEP,s2)
  s2=s2*2.0

  ! d2
  call gaussmf(time,sigma,d20,NSTEP,d2)
  d2=d2*2.0

 do itype=1,size(measurement_type)

 !! initialization 
 adj1(:) = 0.d0 
 adj2(:) = 0.d0

 !! misfit and adjoint source
 call misfit_adj_DD(measurement_type(itype), d1, s1, d2, s2, &
                    NSTEP, deltat, f0, &
                    ntstart1, ntend1, ntstart2, ntend2,&
                    window_type,compute_adjoint, &
                    misfit_output, adj1, adj2)

print*, 'measurement_', trim(measurement_type(itype)) ,' = ', misfit_output
print*, 'squared misfit_', trim(measurement_type(itype)) ,' = ', misfit_output**2
print*

open(IOUT,file=trim(output_dir)//'/adjDD1_'//trim(measurement_type(itype)),status='unknown')
do  i =  ntstart1, ntend1
    write(IOUT,*) time(i), d1(i), s1(i), adj1(i)
enddo
close(IOUT)

open(IOUT,file=trim(output_dir)//'/adjDD2_'//trim(measurement_type(itype)),status='unknown')
do  i =  ntstart2, ntend2
    write(IOUT,*) time(i), d2(i), s2(i), adj2(i)
enddo
close(IOUT)

enddo


end program test_adjoint_sources_DD
