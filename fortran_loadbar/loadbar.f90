!	github.com/itsmoosh/fun/fortran_loadbar/loadbar.f90
!
!	Prints a loading bar to the terminal that is self-overwriting.
!	Please see the README for more information.
!
!	Output cursor is stepped backward after each time the loading
!	bar is updated, unless percent >= 100. Call loadbar(100) after
!	the last sub-100 update in order to place the cursor at a normal
!	location.
!
!	Inputs:
!		percent -- integer in the range [0, 100]
!
!	Author: Marshall 'Moosh' Styczinski
!
!	Last updated: July 18, 2018
!
subroutine loadbar(percent)
	implicit none
	
	integer, intent(in)	:: percent
	integer, parameter	:: mult = 6	! (mult*10)+2 is number of characters used
	character,parameter	:: esc  = char(27)

	integer :: loop		= 1
	integer :: barperc
    integer :: totalbar	= mult*10
	integer :: progress
	integer :: remain

	character*10 :: full10 = '==========', empty10 = '          '
	character*60 :: full = '', empty = ''

	character*6 :: fmtbeg
	character*6 :: fmtend
	character*3 :: progfmt
	character*3 :: remfmt

	barperc = percent

	if(barperc.lt.0) barperc = 0

	if(barperc.ge.100) then
		barperc = 100
		progress = totalbar
		remain = 0
	else
		progress = barperc*totalbar/100
		remain = totalbar - progress - 1
	endif

	do n=1, mult
		full	= trim(full)//full10
		empty	= trim(empty)//empty10
	enddo

	write(progfmt,'(I3)') progress
	fmtbeg   = '(A'//trim(adjustl(progfmt))//'$)'
	write(remfmt,'(I3)') remain
	fmtend = '(A'//trim(adjustl(remfmt))//'$)'

	write(*,'(A1$)') '['
	write(*,fmtbeg) full
	if(progress.lt.totalbar)	write(*,'(A1$)') '>'
	write(*,fmtend) empty
	write(*,'(A1)') ']'
	write(*,'(A,I3,A)') 'Working,', barperc, '% done.'

	if(barperc.lt.100) write(*,'(A,A$)') esc, '[2A'

	return
end subroutine loadbar
