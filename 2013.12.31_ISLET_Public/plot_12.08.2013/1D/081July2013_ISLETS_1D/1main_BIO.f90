!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!     Hogil Kim Memorial Building #501 POSTECH,
!     San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!     Personal site: http://hoangdanhtai.com
!-----------------------------------------------------------------------------------------!
!     08.03.2013: Na=Nb, plot Naa as a funtions of Na=Nb 
!      (thay doi gia tri Nx de so sanh voi kq tinh toan bang tay cua Juyong)
!!    08.03.2013: Chon 1 nguyen tu bat ki (khong phai la NN nhu truoc) de update
!!    14.03.2013: Mo rong cho 3D
!!    18.03.2013: Sua lai phan tinh local fiel, chon random de so sanh ==> OK
!!    20.03.2013: Bo sung SEED
!!    22.03.13: Tinh fluctuation of Paa,Pbb,Pab
!!    27.03.13: Bo sung thay doi Jab
!!    10.05.13: chuyen del_number_a ve integer cho dung
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
      PROGRAM main_bio
      IMPLICIT NONE

      CHARACTER (LEN=150):: CONFIG_INI
      CHARACTER (LEN=50) :: name
      CHARACTER (LEN=15) :: tmp
      CHARACTER (LEN=3)  :: SAt
      CHARACTER (256)    :: Ligne21,Ligne22,Ligne23
     
      REAL    (KIND=8),PARAMETER :: aaa=4. , nul=0.

      INTEGER (KIND=8):: iT,nT,i,j,ip,jp,im,jm,k,km,kp,l
      INTEGER (KIND=8):: natx,naty,natz,natx_p,naty_p,natz_p,number
      INTEGER (KIND=8):: i_loop,i_loop1,i_loop2,n_equi1,n_equi2,n_average,i_times,n_times
      INTEGER (KIND=8) :: inn,jnn,inn_m,inn_p,jnn_m,jnn_p,knn,knn_m,knn_p
      
      INTEGER (KIND=8) :: number_a,na,nb,i_number_a,n_number_a,n_Jab,i_Jab,del_number_a
      
      REAL    (KIND=8) :: T,delT,Tmax,Tmin,x,y,z,n_atom,pro_a_min,pro_a_max
      REAL    (KIND=8) :: rdn_mtp,rdn_configx,rdn_configy,rdn_configz
      REAL    (KIND=8) :: Ha,Hb,Hnna,Hnnb,Jaa,Jbb,Jab,Sa_tmp,Sb_tmp
      REAL    (KIND=8) :: E1_old,E2_old,E1_new,E2_new,energy,E_av
      REAL    (KIND=8) :: Paa,Pbb,Pab,rdn_inn,rdn_jnn,rdn_knn
      REAL    (KIND=8) :: Naa,Nbb,Nab,Naa_av,Nbb_av,Nab_av,Nsa,Nsb,N_total
      REAL    (KIND=8) :: Naa2,Nbb2,Nab2,Naa2_av,Nbb2_av,Nab2_av,FNaa,FNbb,FNab,FPaa,FPbb,FPab
      REAL    (KIND=8) :: Jabmin,Jabmax,del_Jab

      INTEGER (KIND=8),DIMENSION(:),ALLOCATABLE :: tab_ip,tab_im,tab_jp,tab_jm,tab_kp,tab_km

      REAL (KIND=8),DIMENSION(:,:,:),ALLOCATABLE  :: Sa,Sb
                                     
        
!!!=======================================================================================
!!!=======================================================================================
      CALL system('rm -r config_ini_3D')
      CALL system('mkdir config_ini_3D')
      CALL system('rm -r config_3D')
      CALL system('mkdir config_3D')
      CALL system('rm *.dat*')


      CALL ini_rdm_number()
      CALL read_input_parameter()

      !WRITE(*,*)'n_Jaa, Jaamin, Jaamax=',n_Jaa,Jaamin,Jaamax

      natx_p=natx+1 ; naty_p=naty+1 ; natz_p=natz+1

      n_atom=real(natx*naty*natz)

      ALLOCATE(tab_ip(0:natx_p)) ; ALLOCATE(tab_im(0:natx_p))
      ALLOCATE(tab_jp(0:naty_p)) ; ALLOCATE(tab_jm(0:naty_p))
      ALLOCATE(tab_kp(0:natz_p)) ; ALLOCATE(tab_km(0:natz_p))

      ALLOCATE(Sa(-1:natx_p+1,-1:naty_p+1,-1:natz_p+1))
      ALLOCATE(Sb(-1:natx_p+1,-1:naty_p+1,-1:natz_p+1))

      CALL tab_ip_im()

      IF (nT==1) THEN
            delT=0.
      ELSE
            delT=(Tmax-Tmin)/real(nT-1)
      END IF

      IF (n_number_a==1) THEN
            del_number_a=0
      ELSE
            del_number_a=int(n_atom*(pro_a_max-pro_a_min)/real(n_number_a-1))
      END IF

      IF (n_Jab==1) THEN
            del_Jab=0.
      ELSE
            del_Jab=(Jabmax-Jabmin)/real(n_Jab-1)
      END IF

      !WRITE(*,*)'del_Jab=',del_Jab

      WRITE(*,*)"n_number_a",n_number_a,"del_number_a = ",del_number_a
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM ======
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      OPEN(unit=21,file='average_thermal.dat')
      OPEN(unit=22,file='Paa.dat')
      OPEN(unit=23,file='E_itimes.dat')

      DO i_number_a=1,n_number_a
            number_a = int(pro_a_min*n_atom) + del_number_a*(i_number_a-1)
            
            WRITE(*,*)'number_a=',number_a
            
            CALL load_config_ini()
            
            CALL write_config_ini_3D()

            DO i_Jab=1,n_Jab
            Jab=Jabmin+del_Jab*real(i_Jab-1)
                        
            WRITE (*,*)'Jab=',Jab

            DO iT=1,nT
                  WRITE(*,*)'iT = ', iT                        
                  T=Tmin+delT*real(iT-1)
                  
                  CALL equi_lattice1()
                  CALL average_thermal()
       
                  !CALL value_thermal()
                  
                  !WRITE(*,*)'Naa=',Naa,'Nbb=',Nbb,'Nab=',Nab
            
                  CALL write_config_3D()
                  
                  !WRITE(*,*)'E_average=',E_average

            END DO

            END DO

      END DO


      !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      Nsa=0. ; Nsb=0.
      DO k=1,natz
      DO j=1,naty
      DO i=1,natx
            Nsa=Nsa+Sa(i,j,k)
            Nsb=Nsb+Sb(i,j,k)
      END DO      
      END DO
      END DO
      
      WRITE(*,*)'Nsa=',Nsa,'Nsb=',Nsb
            
      !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

      CLOSE(21)    
      CLOSE(22)
      CLOSE(23)

      DEALLOCATE(Sa,Sb)
      DEALLOCATE(tab_ip,tab_im,tab_jp,tab_jm,tab_kp,tab_km)
     
      CONTAINS


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ini_rdm_number()
      IMPLICIT NONE

      INTEGER (KIND=8),DIMENSION(8) :: time
      INTEGER (KIND=8),DIMENSION(50) :: seed

      CALL DATE_AND_TIME(values=time)     ! Get the current time
      seed(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
      CALL RANDOM_SEED(PUT=seed)

      END SUBROUTINE ini_rdm_number 
!!!======================================================================================= 
      SUBROUTINE tab_ip_im()
      IMPLICIT NONE

      DO i=1,natx
            tab_ip(i)=i+1 ; tab_im(i)=i-1
      END DO

      DO j=1,naty
            tab_jp(j)=j+1 ; tab_jm(j)=j-1
      END DO

      DO k=1,natz
            tab_kp(k)=k+1 ; tab_km(k)=k-1
      END DO

      END SUBROUTINE tab_ip_im
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE read_input_parameter() 
!!! OPEN the parameter from file "parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_input_parameter()
      IMPLICIT NONE

      CHARACTER (LEN=150) :: tamp
      OPEN(11,file='1parameter.in')
     
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(A10))')   tamp, CONFIG_INI
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I5))')    tamp, natx
      READ(11, '(A30,(I5))')    tamp, naty
      READ(11, '(A30,(I5))')    tamp, natz
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I5))')    tamp, n_number_a
      READ(11, '(A30,(F12.6))') tamp, pro_a_min
      READ(11, '(A30,(F12.6))') tamp, pro_a_max
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jaa
      READ(11, '(A30,(F12.6))') tamp, Jbb
      READ(11, '(A30,(F12.6))') tamp, Jab
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp, nT
      READ(11, '(A30,(F7.4))')  tamp, Tmin
      READ(11, '(A30,(F7.4))')  tamp, Tmax
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I12))')   tamp,n_equi1
      READ(11, '(A30,(I12))')   tamp,n_equi2
      READ(11, '(A30,(I12))')   tamp,n_average
      READ(11, '(A30,(I12))')   tamp,n_times
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp, n_Jab
      READ(11, '(A30,(F7.4))')  tamp, Jabmin
      READ(11, '(A30,(F7.4))')  tamp, Jabmax
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 

      END SUBROUTINE read_input_parameter


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Initial position configuration
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE load_config_ini()
      IMPLICIT NONE      

      Sa(:,:,:)=0. ; Sb(:,:,:)=0.

!!!------------------------------------------------------
      IF (CONFIG_INI == 'NO') THEN

      na=0 ; nb=0
      
      DO WHILE (na<number_a)
            i=0 ; j=0 ; k=0

            DO WHILE(i==0)
                  CALL random_number(rdn_configx)
                  i=int(rdn_configx*real(natx_p))
            ENDDO

            DO WHILE(j==0)
                  CALL random_number(rdn_configy)
                  j=int(rdn_configy*real(naty_p))
            ENDDO

            DO WHILE(k==0)
                  CALL random_number(rdn_configz)
                  k=int(rdn_configz*real(natz_p))
            ENDDO

            IF (int(Sa(i,j,k))==0) THEN                      
                  Sa(i,j,k)=1.
                  na=na+1
            END IF

      END DO

      DO k=1,natz
      DO j=1,naty
      DO i=1,natx

            IF (int(Sa(i,j,k))==0) THEN   
                  Sb(i,j,k)=1.
                  nb=nb+1
            END IF

      ENDDO
      ENDDO
      ENDDO

      !WRITE(*,*)'na=',na,'nb=',nb

      END IF

!!!------------------------------------------------------
     IF ((CONFIG_INI == 'GS1').or.(CONFIG_INI == 'YES')) THEN

      na=0 ; nb=0

      DO k=1,natz
      DO j=1,naty
      DO i=1,natx

            IF ((int(Sa(i,j,k))==0).and.(na<number_a)) THEN                      
                  Sa(i,j,k)=1.
                  na=na+1
            END IF

      ENDDO
      ENDDO
      ENDDO

      DO k=1,natz
      DO j=1,naty
      DO i=1,natx

            IF (int(Sa(i,j,k))==0) THEN   
                  Sb(i,j,k)=1.
                  nb=nb+1
            END IF

      ENDDO
      ENDDO
      ENDDO

      END IF
      
!!!---------------------------------------------------------------------------------
      IF (CONFIG_INI == 'GS2') THEN

      na=0 ; nb=0


      DO l=1,natx
      
      DO k=1,natz
      DO j=1,naty
      DO i=1,natx

      IF ((int(Sa(i,j,k))==0).and.(na<number_a)) THEN       
            IF ((i==l).or.(i==(natx+1-l)).or.(j==1).or.(j==(naty+1-l))&
            .or.(k==1).or.(k==(natz+1-l))) THEN                      
                  Sa(i,j,k)=1.
                  na=na+1
            END IF       
      END IF
       
      ENDDO
      ENDDO
      ENDDO 
      
      END DO
      
      !!!--------------------
      DO k=1,natz
      DO j=1,naty
      DO i=1,natx

            IF (int(Sa(i,j,k))==0) THEN   
                  Sb(i,j,k)=1.
                  nb=nb+1
            END IF

      ENDDO
      ENDDO
      ENDDO

      END IF

!!!---------------------------------------------------------------------------------
      IF (CONFIG_INI == 'GS3') THEN

      na=0 ; nb=0
      
      DO k=1,natz
      DO j=1,naty
      DO i=1,natx

      IF ((int(Sa(i,j,k))==0).and.(na<number_a)) THEN       
            IF (mod(i+j+k,2)==0) THEN                      
                  Sa(i,j,k)=1.
                  na=na+1
            END IF       
      END IF
       
      ENDDO
      ENDDO
      ENDDO 

      !!!--------------------
      DO k=1,natz
      DO j=1,naty
      DO i=1,natx

            IF (int(Sa(i,j,k))==0) THEN   
                  Sb(i,j,k)=1.
                  nb=nb+1
            END IF

      ENDDO
      ENDDO
      ENDDO

      END IF

      END SUBROUTINE load_config_ini

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_ini_3D()
      IMPLICIT NONE
 
      OPEN(unit=12,file='config_ini_3D/config_ini_3D.pdb')

      DO k=1,natz
      DO j=1,naty
      DO i=1,natx
            
            x=real(i-1)*aaa
            y=real(j-1)*aaa
            z=real(k-1)*aaa
                     
            !IF ((i>5).or.(j>5).or.(k<6)) THEN

            IF (int(Sa(i,j,k))==1) THEN 
                  SAt='Cu'
                  WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x,y,z,nul
            ELSE
                              
            IF (int(Sb(i,j,k))==1) THEN 
                  SAt='H'
                  WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x,y,z,nul

            ELSE
                  SAt='H'
                  WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x,y,z,nul

            END IF
            END IF

            !END IF

      END DO
      END DO
      END DO

      CLOSE(12)

      END SUBROUTINE write_config_ini_3D


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_3D()
      IMPLICIT NONE

      number=10000000+iT*1000+i_number_a+i_Jab-1
      
      WRITE(tmp,'(I8)') number

      name='_BIO_'//TRIM(tmp)
      
      OPEN(unit=13,file='config_3D/config_3D'//trim(name)//'.pdb')

      DO k=1,natz
      DO j=1,naty
      DO i=1,natx

            x=real(i-1)*aaa
            y=real(j-1)*aaa
            z=real(k-1)*aaa
            
            !IF (k==3) THEN
            
            !IF ((i>5).or.(j>5).or.(k<6)) THEN
            
            IF (int(Sa(i,j,k))==1) THEN 
                  SAt='Cu'
                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x,y,z,nul

            ELSE
                              
            IF (int(Sb(i,j,k))==1) THEN 
                  SAt='H'
                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x,y,z,nul

            ELSE
                  SAt='H'
                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x,y,z,nul

            END IF
            END IF

            !END IF

      END DO
      END DO
      END DO

      CLOSE(13)

      END SUBROUTINE write_config_3D


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_lattice1()
      IMPLICIT NONE

      DO i_loop=1,n_equi1
            !CALL equi_lattice()
            CALL value_thermal()
      END DO

      END SUBROUTINE equi_lattice1

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_lattice()
      IMPLICIT NONE
          
      E1_old=0. ; E2_old=0. ; E1_new=0. ; E2_new=0.

      DO k=1,natz
      kp=tab_kp(k)
      km=tab_km(k)
      
      DO j=1,naty
      jp=tab_jp(j)
      jm=tab_jm(j)

      DO i=1,natx
      ip=tab_ip(i)
      im=tab_im(i)
           
           inn=0 ; jnn=0 ; knn=0
           
           DO WHILE (int(Sa(inn,jnn,knn)+Sb(inn,jnn,knn))==0)
           
           CALL random_number(rdn_inn)
           inn=int(rdn_inn*real(natx_p))
           
           CALL random_number(rdn_jnn)
           jnn=int(rdn_jnn*real(naty_p))
           
           CALL random_number(rdn_knn)
           knn=int(rdn_knn*real(natz_p))
                      
           END DO
           CALL value_H()
                
           E1_old=-Sa(i,j,k)*(Jaa*Ha+Jab*Hb)-Sb(i,j,k)*(Jab*Ha+Jbb*Hb)   

!!==================================================================================================           
           IF (int(Sa(i,j,k))/=int(Sa(inn,jnn,knn))) THEN
           
                 CALL value_Hnn()
                 E2_old=-Sa(inn,jnn,knn)*(Jaa*Hnna+Jab*Hnnb)-Sb(inn,jnn,knn)*(Jab*Hnna+Jbb*Hnnb)

                 Sa_tmp=Sa(i,j,k)          ;  Sb_tmp=Sb(i,j,k)
                 Sa(i,j,k)=Sa(inn,jnn,knn) ;  Sb(i,j,k)=Sb(inn,jnn,knn)
                 Sa(inn,jnn,knn)=Sa_tmp    ;  Sb(inn,jnn,knn)=Sb_tmp

                 CALL value_H()
                 CALL value_Hnn()            
                 E1_new=-Sa(i,j,k)*(Jaa*Ha+Jab*Hb)-Sb(i,j,k)*(Jab*Ha+Jbb*Hb) 
                 E2_new=-Sa(inn,jnn,knn)*(Jaa*Hnna+Jab*Hnnb)-Sb(inn,jnn,knn)*(Jab*Hnna+Jbb*Hnnb)

                 CALL random_number(rdn_mtp)

                 IF (exp(-(E1_new+E2_new-E1_old-E2_old)/T) > rdn_mtp) THEN
                 E1_old=E1_new

                 ELSE         
                 Sa(inn,jnn,knn)=Sa(i,j,k) ;  Sb(inn,jnn,knn)=Sb(i,j,k)
                 Sa(i,j,k)=Sa_tmp          ;  Sb(i,j,k)=Sb_tmp

                 END IF

           END IF
!!==================================================================================================

      END DO
      END DO
      END DO
                     
      END SUBROUTINE equi_lattice     
      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_H()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_H()
      IMPLICIT NONE
      
      !!! value H    
      Ha=Sa(im,j,k)+Sa(ip,j,k)+Sa(i,jm,k)+Sa(i,jp,k)+Sa(i,j,kp)+Sa(i,j,km)
      Hb=Sb(im,j,k)+Sb(ip,j,k)+Sb(i,jm,k)+Sb(i,jp,k)+Sb(i,j,kp)+Sb(i,j,km)

      END SUBROUTINE value_H
!!!=======================================================================================
      SUBROUTINE value_Hnn()
      IMPLICIT NONE

      !!! value Hnn     
      inn_m=inn-1 ; inn_p=inn+1 ; jnn_m=jnn-1 ; jnn_p=jnn+1 ; knn_m=knn-1 ; knn_p=knn+1
                    
      Hnna=Sa(inn_m,jnn,knn)+Sa(inn_p,jnn,knn)+Sa(inn,jnn_m,knn)+Sa(inn,jnn_p,knn) &
          +Sa(inn,jnn,knn_m)+Sa(inn,jnn,knn_p)
      
      Hnnb=Sb(inn_m,jnn,knn)+Sb(inn_p,jnn,knn)+Sb(inn,jnn_m,knn)+Sb(inn,jnn_p,knn) &
          +Sb(inn,jnn,knn_m)+Sb(inn,jnn,knn_p)

      END SUBROUTINE value_Hnn

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_thermal()
      IMPLICIT NONE
          
      energy=0. ; E1_old=0. ; E2_old=0. ; E1_new=0. ; E2_new=0.
      Naa=0. ; Nbb=0. ; Nab=0. ; Naa2=0. ; Nbb2=0. ; Nab2=0.

      DO k=1,natz
      kp=tab_kp(k)
      km=tab_km(k)
      
      DO j=1,naty
      jp=tab_jp(j)
      jm=tab_jm(j)

      DO i=1,natx
      ip=tab_ip(i)
      im=tab_im(i)
           
           inn=0 ; jnn=1 ; knn=1
           
           DO WHILE (int(Sa(inn,jnn,knn)+Sb(inn,jnn,knn))==0)
           
           CALL random_number(rdn_inn)
           inn=int(rdn_inn*real(natx_p))
           
           !CALL random_number(rdn_jnn)
           !jnn=int(rdn_jnn*real(naty_p))
           
           !CALL random_number(rdn_knn)
           !knn=int(rdn_knn*real(natz_p))
                      
           END DO

           CALL value_H()
                
           E1_old=-Sa(i,j,k)*(Jaa*Ha+Jab*Hb)-Sb(i,j,k)*(Jab*Ha+Jbb*Hb)   

!!========================================================================================        
            IF (int(Sa(i,j,k))/=int(Sa(inn,jnn,knn))) THEN           

                 CALL value_Hnn()
                 E2_old=-Sa(inn,jnn,knn)*(Jaa*Hnna+Jab*Hnnb)-Sb(inn,jnn,knn)*(Jab*Hnna+Jbb*Hnnb)

                 Sa_tmp=Sa(i,j,k)          ;  Sb_tmp=Sb(i,j,k)
                 Sa(i,j,k)=Sa(inn,jnn,knn) ;  Sb(i,j,k)=Sb(inn,jnn,knn)
                 Sa(inn,jnn,knn)=Sa_tmp    ;  Sb(inn,jnn,knn)=Sb_tmp

                 CALL value_H()
                 CALL value_Hnn()            
                 E1_new=-Sa(i,j,k)*(Jaa*Ha+Jab*Hb)-Sb(i,j,k)*(Jab*Ha+Jbb*Hb) 
                 E2_new=-Sa(inn,jnn,knn)*(Jaa*Hnna+Jab*Hnnb)-Sb(inn,jnn,knn)*(Jab*Hnna+Jbb*Hnnb)

                 CALL random_number(rdn_mtp)

                 IF (exp(-(E1_new+E2_new-E1_old-E2_old)/T) > rdn_mtp) THEN
                 E1_old=E1_new

                 ELSE         
                 Sa(inn,jnn,knn)=Sa(i,j,k) ;  Sb(inn,jnn,knn)=Sb(i,j,k)
                 Sa(i,j,k)=Sa_tmp          ;  Sb(i,j,k)=Sb_tmp

                 END IF

           END IF
!!==========================================================================================                       
           
           energy=energy+E1_old

            !!! Tinh nAA, NAB,...-----------------------------------
            CALL value_H()
            Naa=Naa+Sa(i,j,k)*Ha ; Nbb=Nbb+Sb(i,j,k)*Hb ; Nab=Nab+Sa(i,j,k)*Hb+Sb(i,j,k)*Ha
            !!!-----------------------------------------------------

      END DO
      END DO
      END DO
      

      energy=energy/2./n_atom

      Naa=Naa/2. ; Nbb=Nbb/2. ; Nab=Nab/2.
      
      Naa2=Naa**2.; Nbb2=Nbb**2. ; Nab2=Nab**2.
      
      END SUBROUTINE value_thermal


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_Naa()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_Naa()
      IMPLICIT NONE
          
      Naa=0. ; Nbb=0. ; Nab=0.

      DO k=1,natz
      kp=tab_kp(k)
      km=tab_km(k)

      DO j=1,naty
      jp=tab_jp(j)
      jm=tab_jm(j)    
      DO i=1,natx
      ip=tab_ip(i)
      im=tab_im(i)
            
            CALL value_H()
                   
            Naa=Naa+Sa(i,j,k)*Ha ; Nbb=Nbb+Sb(i,j,k)*Hb ; Nab=Nab+Sa(i,j,k)*Hb+Sb(i,j,k)*Ha

      END DO
      END DO
      END DO
      
      Naa=Naa/2. ; Nbb=Nbb/2. ; Nab=Nab/2.

      END SUBROUTINE value_Naa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE average_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE average_thermal()

      IMPLICIT NONE

      E_av=0.; Naa_av=0.; Nbb_av=0. ; Nab_av=0.
      Naa2_av=0.; Nbb2_av=0. ; Nab2_av=0.

      DO i_times=1,n_times

            DO i_loop1=1,n_equi2
                  CALL equi_lattice()
            END DO

            DO i_loop2=1,n_average
                  CALL value_thermal()        
                  
                  E_av=E_av+energy
                  Naa_av=Naa_av+Naa ; Nbb_av=Nbb_av+Nbb ; Nab_av=Nab_av+Nab
                  Naa2_av=Naa2_av+Naa2 ; Nbb2_av=Nbb2_av+Nbb2 ; Nab2_av=Nab2_av+Nab2
                                      
            END DO
      
      WRITE(Ligne23,*) i_times,E_av/real(i_times*n_average)
      WRITE(23,'(a)') trim(Ligne23)

      END DO

      E_av=E_av/real(n_times*n_average)
      
      Naa_av=Naa_av/real(n_times*n_average)
      Nbb_av=Nbb_av/real(n_times*n_average)
      Nab_av=Nab_av/real(n_times*n_average)
      
      Naa2_av=Naa2_av/real(n_times*n_average)
      Nbb2_av=Nbb2_av/real(n_times*n_average)
      Nab2_av=Nab2_av/real(n_times*n_average)      
       
      FNaa=Naa2_av-Naa_av**2.
      FNbb=Nbb2_av-Nbb_av**2.
      FNab=Nab2_av-Nab_av**2.      
      
      N_total=(Naa_av+Nbb_av+Nab_av)
            
      Paa=Naa_av/N_total ; Pbb=Nbb_av/N_total ; Pab=Nab_av/N_total
      
      FPaa=FNaa/N_total ; FPbb=FNbb/N_total ; FPab=FNab/N_total

      WRITE(Ligne21,*) Jab,(1.-real(number_a/n_atom)),E_av,Naa_av,Nbb_av,Nab_av,FNaa,FNbb,FNab
      WRITE(21,'(a)') trim(Ligne21)
      
      WRITE(Ligne22,*) T,(1.-real(number_a/n_atom)),Jab,Paa,Pbb,Pab,FPaa,FPbb,FPab
      WRITE(22,'(a)') trim(Ligne22)
      
      END SUBROUTINE average_thermal

      END PROGRAM main_bio
      

      
