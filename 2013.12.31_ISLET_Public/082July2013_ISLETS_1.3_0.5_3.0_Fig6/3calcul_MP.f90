!!! ======================================================================================
!!! ====================================================================================== 
!!! Chuong trinh tong hop ket qua tu cach tinh MP
!!! 22.03.2013: 
!!! ======================================================================================
!!! ======================================================================================      
      PROGRAM calcul_MP
      IMPLICIT NONE
      
      INTEGER (KIND=4), PARAMETER :: nF=16, n_number_a=21

      INTEGER (KIND=4) :: i,j
    
      REAL    (KIND=8),DIMENSION(n_number_a) :: T,Pa,na,Paa_MP,Pbb_MP,Pab_MP           
      REAL    (KIND=8),DIMENSION(nF,n_number_a) :: Paa,Pbb,Pab
     
!!!!######################################################################################
!!!!######################################################################################
!!!!######################################################################################

      
!!!---------------------------------------------------------------------------------------     
!!! Doc gia tri Paa,Pbb,Pab tu cac file       
      DO i=1,nF
            CALL read_data()
    
            DO j=1,n_number_a                 
                  READ(12,*)T(j),Pa(j),na(j),Paa(i,j),Pbb(i,j),Pab(i,j)
            END DO
            
      END DO 

      CLOSE(12) 
                 
!!!---------------------------------------------------------------------------------------
!!! Tinh gia tri trung binh Paa,Pbb,Pab

      OPEN(unit=22,file='Paa_MP.dat') 
       
      Paa_MP(:)=0. ; Pbb_MP(:)=0. ; Pab_MP(:)=0.
      
      DO j=1,n_number_a

            DO i=1,nF
                  Paa_MP(j)=Paa_MP(j)+Paa(i,j)
                  Pbb_MP(j)=Pbb_MP(j)+Pbb(i,j)
                  Pab_MP(j)=Pab_MP(j)+Pab(i,j)
                     
            END DO
            
            Paa_MP(j)=Paa_MP(j)/real(nF)
            Pbb_MP(j)=Pbb_MP(j)/real(nF)
            Pab_MP(j)=Pab_MP(j)/real(nF)

            WRITE(22,*)T(j),Pa(j),na(j),Paa_MP(j),Pbb_MP(j),Pab_MP(j)

      END DO
     

      CLOSE(22)


      CONTAINS
      
!!!!######################################################################################  
!!!!######################################################################################
!!!!######################################################################################

     
!!!=======================================================================================
!!!=======================================================================================
!!!  SUBROUTINE read_data()
!!!=======================================================================================

      SUBROUTINE read_data()
      IMPLICIT NONE
     
      IF (i==1) THEN
            OPEN(unit=12,file='1/Paa.dat')
      END IF   

      IF (i==2) THEN
            OPEN(unit=12,file='2/Paa.dat')
      END IF 

      IF (i==3) THEN
            OPEN(unit=12,file='3/Paa.dat')
      END IF

      IF (i==4) THEN
            OPEN(unit=12,file='4/Paa.dat')
      END IF 

      IF (i==5) THEN
            OPEN(unit=12,file='5/Paa.dat')
      END IF

      IF (i==6) THEN
            OPEN(unit=12,file='6/Paa.dat')
      END IF

      IF (i==7) THEN
            OPEN(unit=12,file='7/Paa.dat')
      END IF 
      
      IF (i==8) THEN
            OPEN(unit=12,file='8/Paa.dat')
      END IF

      IF (i==9) THEN
            OPEN(unit=12,file='9/Paa.dat')
      END IF

      IF (i==10) THEN
            OPEN(unit=12,file='10/Paa.dat')
      END IF

      IF (i==11) THEN
            OPEN(unit=12,file='11/Paa.dat')
      END IF   

      IF (i==12) THEN
            OPEN(unit=12,file='12/Paa.dat')
      END IF 

      IF (i==13) THEN
            OPEN(unit=12,file='13/Paa.dat')
      END IF

      IF (i==14) THEN
            OPEN(unit=12,file='14/Paa.dat')
      END IF 

      IF (i==15) THEN
            OPEN(unit=12,file='15/Paa.dat')
      END IF

      IF (i==16) THEN
            OPEN(unit=12,file='16/Paa.dat')
      END IF

      IF (i==17) THEN
            OPEN(unit=12,file='17/Paa.dat')
      END IF 
      
      IF (i==18) THEN
            OPEN(unit=12,file='18/Paa.dat')
      END IF

      IF (i==19) THEN
            OPEN(unit=12,file='19/Paa.dat')
      END IF

      IF (i==20) THEN
            OPEN(unit=12,file='20/Paa.dat')
      END IF

      
      END SUBROUTINE read_data
!!!=======================================================================================
     
      END PROGRAM calcul_MP


