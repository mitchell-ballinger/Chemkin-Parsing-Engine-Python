program file_op
implicit none
	integer, parameter :: m = 2**16, n = 512, q = 64
	integer :: i,j,C, LENGTH, LENGTH_A, LENGTH_B, LENGTH_C, EOF, nrct, ndis, nsp, INTEG
	integer :: INT_BUFFER, FLAG
	real :: REAL_BUFFER
	character(16) :: STR_BUFFER, FOLDER
	character(n), dimension(m) :: LINE
	integer, dimension(m,q) :: SPLITS
	real, dimension(q,m) :: VALUES, tbe
	real, dimension(q) :: wt
	real, dimension(6,m) :: nd, stc
	real, dimension(m,4) :: ctroe
	real, dimension(3,m) :: rctdata, rctdatatroe
	CHARACTER(100) :: ARG
	
!=============================================================================
!==== COMMAND LINE INTERFACE
!=============================================================================
	CALL GETARG(1,ARG)
	if (trim(ARG) == '' .or. trim(ARG) == '-h') then
		WRITE(*,*) 'Program Use: read.exe -f <FOLDERNAME>'
		stop
	end if
	if (trim(ARG) == '-f' .or. trim(ARG) == '-F') then
		CALL GETARG(2,ARG)
		FOLDER = trim(ARG)
	end if

!=============================================================================
!==== Read CONFIG File : <NUM_RXNS, NUM_SPECIES, NUM_SCALES>
!=============================================================================
	open (2, file= (trim(FOLDER) // '/config.csv'))
	read(2,*) nrct, nsp
	close(2)
	
!=============================================================================
!==== Read WEIGHTS File : <NUM_RXNS, NUM_SPECIES, NUM_SCALES>
!=============================================================================
	open (2, file= (trim(FOLDER) // '/weights.csv'))
	do i=1,nsp
		read(2,*) wt(i)
	end do
	close(2)
	
!=========================================================================================
!==== Read Three Body File : <[RXN]:6, [STO]:6, [RATE]:3, [TROE]:3, [SCALES]: ndis>
!==== Read Hybrid File : <[RXN]:6, [STO]:6, [LOW]:3, [RATE]:3, [TROE]:3, [SCALES]: ndis>
!==== Read Exchange File : <[RXN]:6, [STO]:6, [RATE]:3>
!=========================================================================================
	open (3, file= (trim(FOLDER) // '/three_body.csv'))
	LENGTH_A = 0
	LENGTH_B = 0
	LENGTH_C = 0
	do i=1,nrct
		read(3,'(A)',IOSTAT=EOF) LINE(i)
		IF (EOF < 0) EXIT
		LENGTH_A = i
	end do
	close(3)
	open (3, file= (trim(FOLDER) // '/hybrid.csv'))
	do i=1,nrct
		read(3,'(A)',IOSTAT=EOF) LINE(i+LENGTH_A)
		IF (EOF < 0) EXIT
		LENGTH_B = i
	end do
	close(3)
	open (3, file= (trim(FOLDER) // '/exchange.csv'))
	do i=1,nrct
		read(3,'(A)',IOSTAT=EOF) LINE(i+LENGTH_A+LENGTH_B)
		IF (EOF < 0) EXIT
		LENGTH_C = i
	end do
	close(3)
	
	ndis = LENGTH_A + LENGTH_B ! Number of M-Type Reactions with Scales
	
	do i=1,nrct
	SPLITS(i,1) = 0
	C = 2
	do j=1,len(LINE(i))
		if (LINE(i)(j:j) == ' ') then 
			SPLITS(i,C) = j
			exit
		end if
		if (LINE(i)(j:j) == ',') then
			SPLITS(i,C) = j
			C = C + 1
		end if
	end do
	do j=1,q
		if (SPLITS(i,j+1) == 0) exit
		STR_BUFFER = LINE(i)(SPLITS(i,j)+1:SPLITS(i,j+1)-1)
		read(STR_BUFFER,'(EN20.0E2)') REAL_BUFFER
		VALUES(j,i) =  REAL_BUFFER
	end do
	end do
	
	do i=1,LENGTH_A !Three Body
		nd(1:6,i) = VALUES(1:6,i)
		stc(1:6,i) = VALUES(7:12,i)
		rctdata(1:3,i) = VALUES(13:15,i)
		ctroe(i,1:4) = VALUES(16:19,i)
		tbe(1:nsp,i) = VALUES(20:20+nsp,i)
	end do

	do i=(LENGTH_A+1),(LENGTH_A+LENGTH_B) !Hybrid
		nd(1:6,i) = VALUES(1:6,i)
		stc(1:6,i) = VALUES(7:12,i)
		rctdata(1:3,i) = VALUES(13:15,i)
		rctdatatroe(1:3,i) = VALUES(16:18,i)
		ctroe(i,1:4) = VALUES(19:22,i)
		tbe(1:nsp,i) = VALUES(23:23+nsp,i) 
	end do

	do i=(LENGTH_A+LENGTH_B+1),nrct !Exchange
		nd(1:6,i) = VALUES(1:6,i)
		stc(1:6,i) = VALUES(7:12,i)
		rctdata(1:3,i) = VALUES(13:15,i) 
	end do
!====================================================================
!===== Write Scales Full
!====================================================================
	do i=1,nsp
		write(*,*) tbe(i,1:(LENGTH_A+LENGTH_B))
	end do
!====================================================================
!===== Write Scales EXC
!====================================================================
!	do i=1,nsp
!	FLAG = 0
!	do j=1,ndis
!		if (tbe(i,j) == 1.0) FLAG = FLAG + 1
!	end do
!		if (FLAG .ne. ndis) write(*,*) tbe(i,1:(LENGTH_A+LENGTH_B))
!	end do
end 