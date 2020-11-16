PROGRAM dipolef
!
! convert dipole matrix file from binary to text
!
   IMPLICIT NONE
   
   INTEGER, PARAMETER :: dp = selected_real_kind(14,200)
   INTEGER :: nbnd, nks, nspin, ib1, ib2, ik
   INTEGER :: iun = 28
   
   COMPLEX(dp), ALLOCATABLE :: mat(:,:,:,:)
   
   OPEN(iun, FILE='edipole', FORM='UNFORMATTED')
   READ(iun) nbnd, nks, nspin
   
   ALLOCATE( mat(nbnd, nbnd, nks, 3) )
   READ(iun) mat
   CLOSE(iun)
   
   OPEN(iun, FILE='edipolef', FORM='FORMATTED')
   WRITE(iun, '(3I6)') nspin, nks, nbnd
   
   
      DO ik = 1, nks
         DO ib1 = 1, nbnd
            DO ib2 = 1, nbnd
               WRITE(iun,'(2I6,3(2F20.10))') ib2, ib1, &
                      mat(ib2, ib1, ik, 1:3)
            ENDDO
         ENDDO
      ENDDO
      
   CLOSE(iun)
      
END PROGRAM
