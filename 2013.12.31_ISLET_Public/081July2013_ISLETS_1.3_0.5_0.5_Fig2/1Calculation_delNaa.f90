      PROGRAM cal_delpaa
      IMPLICIT NONE
      
      CHARACTER (256)    :: Ligne21
      
      INTEGER (KIND=8),PARAMETER :: n=51
      
      INTEGER (KIND=8) :: i
      
      REAL (KIND=8) :: E,Pb
      
      
      REAL (KIND=8),DIMENSION(:),ALLOCATABLE  :: Jab,Naa,Nbb,Nab,delNaa,delNbb,delNab
      
      ALLOCATE(Jab(0:n),Naa(0:n),Nbb(0:n),Nab(0:n),delNaa(0:n),delNbb(0:n),delNab(0:n))
            
      
!!!!======================================================================================
!!!!======================================================================================
      CALL system('rm delNaa.dat*')  
      OPEN(unit=12,file='a.dat')
      OPEN(unit=21,file='delNaa.dat')

      Naa(:)=0. ; Nbb(:)=0. ; Nab(:)=0.
      
      DO i=1,n
            !READ(12,*)T,Pb,Jab(i),Paa(i),Pbb(i),Pab(i)

            READ(12,*)Jab(i),Pb,E,Naa(i),Nbb(i),Nab(i)            

            delNaa(i)=Naa(i)-Naa(i-1)
            delNbb(i)=Nbb(i)-Nbb(i-1)
            delNab(i)=Nab(i)-Nab(i-1)
                  
      END DO
      
      DO i=2,n
           
            WRITE(Ligne21,*)Jab(i),-delNaa(i),-delNbb(i),delNab(i),-delNbb(i)/delNab(i)
            WRITE(21,'(a)') trim(Ligne21)
                  
      END DO
      
     !!!WRITE(Ligne21,*) Jab,(1.-real(number_a/n_atom)),E_av,Naa_av,Nbb_av,Nab_av,FNaa,FNbb,FNab
     !!! WRITE(Ligne22,*) T,(1.-real(number_a/n_atom)),Jab,Paa,Pbb,Pab,FPaa,FPbb,FPab
      
      CLOSE(12)
      CLOSE(21)    
    
      END PROGRAM cal_delpaa          
