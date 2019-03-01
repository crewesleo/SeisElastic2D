module myString_mod
implicit none

INTERFACE VAL2A     ! value to ascii
  module procedure i2a, f2a
END INTERFACE
INTERFACE filename  ! construct file name
  module procedure filename0, filename1, filename11, filename2, filename3
END INTERFACE
contains
!==========================

function i2a(num)   ! integer to ascii
integer, intent(in) :: num
character(len=20) i2a
write(i2a,'(I20)')  num
i2a = adjustl(i2a)   ! left align, removing leading spaces
end function
!---------------
function f2a(num)   ! float to ascii
real, intent(in)    :: num
character(len=20) f2a
write(f2a,'(F20.3)')  num
f2a = adjustl(f2a)
end function
!***************
SUBROUTINE filename0(fname,str1,num)
character(len=*), intent(in) :: str1
integer, intent(in) :: num
character(len=*),intent(out):: fname  !!!
! construct a filename by concatebating: str1, atoi(num)
fname = trim(str1) // trim(val2a(num))
END SUBROUTINE
!---------------
SUBROUTINE filename1(fname,str1,num,str2)
character(len=*), intent(in) :: str1, str2
integer, intent(in) :: num
character(len=*),intent(out):: fname  !!!
! construct a filename by concatebating: str1, atoi(num), str2
fname = trim(str1) // trim(val2a(num)) // trim(str2)
END SUBROUTINE
!---------------
SUBROUTINE filename11(fname,str1,num,num2)
character(len=*), intent(in) :: str1
integer, intent(in) :: num, num2
character(len=*),intent(out):: fname  !!!
fname = trim(str1) // trim(val2a(num)) // 'x' // trim(val2a(num2))
END SUBROUTINE
!---------------
SUBROUTINE filename2(fname,str1,num1,num2,str2)
character(len=*), intent(in) :: str1, str2
integer, intent(in) :: num1, num2
character(len=*),intent(out):: fname !!!
! construct a filename by concatebating: 
!   str1, atoi(num1), 'x', atoi(num2), str2
fname = trim(str1) // trim(val2a(num1)) // 'x' // trim(val2a(num2)) //trim(str2)
END SUBROUTINE
!---------------
SUBROUTINE filename3(fname,str1,num1,str2,num2,str3)
character(len=*), intent(in) :: str1, str2, str3
integer, intent(in) :: num1, num2
character(len=*),intent(out):: fname !!!
fname = trim(str1)//trim(val2a(num1))//'x'//trim(str2)//trim(val2a(num2))//trim(str3)
END SUBROUTINE
!=========================
SUBROUTINE fileparts(pathstr, fname, ext, fullname)
! like the MATLAB function w/ the same name
! e.g., fullname = 'H:\user4\matlab\classpath.txt';  % produces ==>
! pathstr = H:\user4\matlab
! fname   = classpath
! ext     = .txt
!    Windows: '\';  UNIX: '/'
!    Without this separator, the name is taken as a file name
! e.g., fullname = '\.cshrc';  % ==>
! pathstr = \   % important, because root dir has a special meaning.  
!               % If leave it as empty, this would mean current directory
! fname   =
! ext     = .cshrc

character(len=*), intent(out):: pathstr, fname, ext
character(len=*), intent(in) :: fullname
integer   ::  loc1, loc2
loc1  = scan(fullname, '\', BACK = .TRUE.)
loc2  = scan(fullname, '/', BACK = .TRUE.)
if (loc1+loc2 == 0) then ! both 0 ==> no dir separator
  pathstr= ''
  fname  = fullname
else ! one of them is non-zero
  if (loc1==0) then
    loc1 = loc2 ! always make loc1 the desired location
  else if (loc1*loc2 /=0) then ! both non-zero
    loc1 = max(loc1, loc2)
  end if ! done bookkeeping
  
  if (loc1==1) then
    pathstr = '/'
  else
    pathstr = fullname(1:loc1-1)
  end if
  fname = fullname(loc1+1:)
end if  ! done pathstr
  
loc1 = scan(fname, '.', BACK = .TRUE.)
if (loc1==0) then
  ext = ''
else
  ext = fname(loc1:)
  fname = fname(:loc1-1)
end if

END SUBROUTINE
!--------------
SUBROUTINE dir_insert(fullname, num)
! fname = dir\abc.bin; num = 4
! fname = dir4\abc.bin
character(len=*), intent(inout):: fullname
integer,          intent(in)   :: num
character(len=150) :: dir
character(len=50)  :: filenm
character(len=8)   :: ext
call fileparts(dir, filenm, ext, fullname)
fullname = trim(dir)// trim(val2a(num))// '/' // trim(filenm) // trim(ext)
END SUBROUTINE
!--------------
SUBROUTINE filename_insert(fname,str1,num1)
! E.g.,
! fname: abc.bin  (ext must be shorter than 4 char's)
! str:   NW
! num1:  4  ==>
! fname: abcNW4.bin
character(len=*), intent(inout):: fname
character(len=*), intent(in)   :: str1
integer,          intent(in)   :: num1
integer   :: loc
loc  = scan(fname, '.', BACK = .TRUE.)
if ( loc==0 .or. (LEN_TRIM(fname)-loc > 3) ) then
  fname = trim(fname)//trim(str1)//trim(val2a(num1))
else
  fname = fname(1:loc-1)//trim(str1)//trim(val2a(num1))//'.'//trim(fname(loc+1:))
end if
END SUBROUTINE

!==========================
end module

