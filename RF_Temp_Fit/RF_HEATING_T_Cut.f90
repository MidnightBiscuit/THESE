!~ To improve speed:
!~ - on the two simulations where dt is cte, then use compilation directives to make them constants!
!~ - make constant the ia ib indexes for the distribution of the ions across the cores
!~ - make private rx,rz,ry on the coulomb calculation
!~ - merge the accelartion and the velocity update openmp (to be tested if it helps, in the KZM code it did far worst)

! code taken from 20210415
! I adapt it to save all info at each simul time step

include 'mkl_vsl.f90'

program main
USE MKL_VSL_TYPE
USE MKL_VSL
use omp_lib
implicit none

integer :: flag_exit
!1 = Save the ion cloud positions (nothing else) at the same time than the Temperature ; 0 = Do not save it.
!dir$ define save_trj_cloud = 0
!dir$ if (save_trj_cloud .eq. 1)
real, allocatable, dimension(:,:,:) ::  save_trj_Cloud
!dec$ endif

include 'variable_definition.f90'

double precision, parameter :: T_threshold = 100* 3/m_kb_x_inv_n_ions! Kelvin 

double precision, parameter    :: cte_dE_dt(2) = (/ charge*Omega*V_rf/r_0**2, mass*wz2*0.5d+00  /)

    time00 = omp_get_wtime()
    call distribute_ions()

	print*, ''
    print*, '*****************************************************************'
    print*, ' ==>  Welcome in RF_relax_with_Langevin_init_Quench24_2.f90  <== '
    print*, '            Generating initial Ion Cloud from file               '
    print*, '*****************************************************************'

    call prepare_cold_ion_cloud_from_file()

    print*, 'Total time', omp_get_wtime() - time00
contains


subroutine prepare_cold_ion_cloud_from_file()
implicit none
    print*, 'dt', dt
!~     pause
    call create_files_names()
    call init_from_file_Langevin()
!~     call manual_init_1ion()
    call initialize_cos_sin()

    j_end   = (j_start + i_fly_Verlet*n_dt) - modulo(j_start + i_fly_Verlet*n_dt,n_dt) + int(0.25*n_dt)
    print*, 'j_end', j_end
    allocate(save_temperature(floor((j_end-j_start) / real(n_save_temp)) + 10, 0:8))
	!dir$ if(save_trj_cloud .eq. 1)
		allocate( save_trj_Cloud(n_ions,6, floor((j_end-j_start) / real(n_save_temp)) + 10))
	!dir$ endif

!~     size_Energy_array = floor((j_end-j_start) / real(n_dt)) + 10
    size_Energy_array = j_end-j_start + 1
    allocate(save_Energies(size_Energy_array,3))

    call fly()
    call save_data()
!~     call save_Energies_to_file()
endsubroutine

subroutine fly()
implicit none
    print*, ''
    print*, 'Start Laser Fly'
    Ep = 0.0d00; Ec = 0.0d00; Ek = 0.0d00; j_save_E = 0
    vdt = dt; vdt1 = dt1; vdt2 = dt2;
    flag_exit = 0
    do jj = j_start , j_end
        call update_step_fix_dt_measuring_Energies()
        call Measure_of_Temperature()
        call Check_save_temp2()
        if (flag_exit == 1) then
			print*, 'T Threshold reached !'
			exit
		endif
    enddo
    j_start = jj;
    print*, 'Finished Cool Fly'
endsubroutine

subroutine update_step_fix_dt()
implicit none
double precision, dimension(ni):: Ek_aux, Ec_aux, Ux_aux, Uy_aux, Ec_aux2,dE_dt_aux
double precision, dimension(3) :: rji
double precision, dimension(2) :: cte_aux
double precision, parameter    :: a_cte_P(4) = (/charge*V_st / r_0**2  - 0.25*mass*wz2, -charge*V_st / r_0**2  - 0.25*mass*wz2, charge*V_rf / r_0**2, 0.5*mass*wz2 /)
double precision               :: cos_jj, sin_jj,r2inv, Ux, Uy,Uz
integer                        :: i, im, j

    cos_jj = cosi(mod(jj,n_dt)+1)
    cte_aux(1) = a_cte_Quad_RF_LC(1) + a_cte_Quad_RF_LC(3)*cos_jj
    cte_aux(2) = a_cte_Quad_RF_LC(2) - a_cte_Quad_RF_LC(3)*cos_jj
!-----------------------------------------------------------------------
    r_x = r_x + v_x*dt + a1_x*dt2;
    r_y = r_y + v_y*dt + a1_y*dt2;
    r_z = r_z + v_z*dt + a1_z*dt2;

    !$omp parallel default(none) private(im, r2inv, rji,i,j) shared (a2_x,a2_y,a2_z, Ec_aux) firstprivate(ia,ib,cte_aux,r_x,r_y,r_z)
        im = omp_get_thread_num()+1
        a2_x(ia(im):ib(im)) =          cte_aux(1) * r_x(ia(im):ib(im));
        a2_y(ia(im):ib(im)) =          cte_aux(2) * r_y(ia(im):ib(im));
        a2_z(ia(im):ib(im)) = a_cte_Quad_RF_LC(4) * r_z(ia(im):ib(im));
        do i = ia(im), ib(im)
            do j = 1, n_ions
                rji(1)  = r_x(j) - r_x(i)
                rji(2)  = r_y(j) - r_y(i)
                rji(3)  = r_z(j) - r_z(i)
                r2inv   = 1.d0/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3)+softening)
                r2inv   = r2inv * r2inv * r2inv * alpha
                a2_x(i) = a2_x(i) - rji(1)*r2inv
                a2_y(i) = a2_y(i) - rji(2)*r2inv
                a2_z(i) = a2_z(i) - rji(3)*r2inv
            enddo
        enddo
    !$omp end parallel

    v_x = v_x + dt1*(a1_x + a2_x); a1_x = a2_x;
    v_y = v_y + dt1*(a1_y + a2_y); a1_y = a2_y;
    v_z = v_z + dt1*(a1_z + a2_z); a1_z = a2_z;

    dE_dt = cte_dE_dt(1)*sini(mod(jj,n_dt)+1)*(sum(r_x*r_x) - sum(r_y*r_y)) - cte_dE_dt(2)*sum(r_x*v_x + r_y*v_y - 2*r_z*v_z)
    
    j_save_E = j_save_E + 1
    save_Energies(j_save_E,1) = Ep
    save_Energies(j_save_E,2) = Ec
    save_Energies(j_save_E,3) = Ek
    save_Energies(j_save_E,4) = dE_dt
    save_Energies(j_save_E,5) = mass/(kb*real(n_ions)) *sum(v_x(:n_ions)**2)  ! m_kb_x_inv_n_ions *save_temperature(i,4:6)
    save_Energies(j_save_E,6) = mass/(kb*real(n_ions)) *sum(v_y(:n_ions)**2)
    save_Energies(j_save_E,7) = mass/(kb*real(n_ions)) *sum(v_z(:n_ions)**2)
end subroutine

subroutine update_step_fix_dt_measuring_Energies()
implicit none
double precision, dimension(ni):: Ek_aux, Ec_aux, Ux_aux, Uy_aux, Ec_aux2,dE_dt_aux
double precision, dimension(3) :: rji
double precision, dimension(2) :: cte_aux
double precision, parameter    :: a_cte_P(4) = (/charge*V_st / r_0**2  - 0.25*mass*wz2, -charge*V_st / r_0**2  - 0.25*mass*wz2, charge*V_rf / r_0**2, 0.5*mass*wz2 /)
double precision               :: cos_jj, sin_jj,r2inv, Ux, Uy,Uz
integer                        :: i, im, j

    cos_jj = cosi(mod(jj,n_dt)+1)
    cte_aux(1) = a_cte_Quad_RF_LC(1) + a_cte_Quad_RF_LC(3)*cos_jj
    cte_aux(2) = a_cte_Quad_RF_LC(2) - a_cte_Quad_RF_LC(3)*cos_jj
!-----------------------------------------------------------------------
    r_x = r_x + v_x*dt + a1_x*dt2;
    Ux = sum(r_x*r_x)

    r_y = r_y + v_y*dt + a1_y*dt2;
    Uy = sum(r_y*r_y)

    r_z = r_z + v_z*dt + a1_z*dt2;
    Uz = sum(r_z*r_z)

    cos_jj = a_cte_P(3)*cos_jj
    Ep = (a_cte_P(1) - cos_jj)*Ux + (a_cte_P(2) + cos_jj)*Uy + a_cte_P(4)*Uz;

    !$omp parallel default(none) private(im, r2inv, rji,i,j) shared (a2_x,a2_y,a2_z, Ec_aux) firstprivate(ia,ib,cte_aux,r_x,r_y,r_z)
        im = omp_get_thread_num()+1
        a2_x(ia(im):ib(im)) =          cte_aux(1) * r_x(ia(im):ib(im));
        a2_y(ia(im):ib(im)) =          cte_aux(2) * r_y(ia(im):ib(im));
        a2_z(ia(im):ib(im)) = a_cte_Quad_RF_LC(4) * r_z(ia(im):ib(im));
        Ec_aux(im) = 0.0d+00
        do i = ia(im), ib(im)
            do j = 1, n_ions
            if (i.ne.j) then
                rji(1)  = r_x(j) - r_x(i)
                rji(2)  = r_y(j) - r_y(i)
                rji(3)  = r_z(j) - r_z(i)
                r2inv   = 1.d0/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3))
                Ec_aux(im) = Ec_aux(im) + r2inv
                r2inv   = r2inv * r2inv * r2inv * alpha
                a2_x(i) = a2_x(i) - rji(1)*r2inv
                a2_y(i) = a2_y(i) - rji(2)*r2inv
                a2_z(i) = a2_z(i) - rji(3)*r2inv
            endif
            enddo
        enddo
    !$omp end parallel
    Ec = sum(Ec_aux)

    v_x = v_x + dt1*(a1_x + a2_x); a1_x = a2_x;
    v_y = v_y + dt1*(a1_y + a2_y); a1_y = a2_y;
    v_z = v_z + dt1*(a1_z + a2_z); a1_z = a2_z;

    Ek = sum( v_x*v_x + v_y*v_y + v_z*v_z )
    S_rf_x  = cte_dE_dt(1)*sini(mod(jj,n_dt)+1)*(Ux)
    S_rf_y  = -cte_dE_dt(1)*sini(mod(jj,n_dt)+1)*(Uy)
    dE_dt = cte_dE_dt(1)*sini(mod(jj,n_dt)+1)*(Ux - Uy) - cte_dE_dt(2)*sum(r_x*v_x + r_y*v_y)
    
    j_save_E = j_save_E + 1
    ! save_Energies(j_save_E,1) = Ep
    ! save_Energies(j_save_E,2) = Ec
    ! save_Energies(j_save_E,3) = Ek
    ! save_Energies(j_save_E,4) = dE_dt
    ! save_Energies(j_save_E,5) = mass/(kb*real(n_ions)) *sum(v_x(:n_ions)**2)  ! m_kb_x_inv_n_ions *save_temperature(i,4:6)
    ! save_Energies(j_save_E,6) = mass/(kb*real(n_ions)) *sum(v_y(:n_ions)**2)
    ! save_Energies(j_save_E,7) = mass/(kb*real(n_ions)) *sum(v_z(:n_ions)**2)
	
	save_Energies(j_save_E,1) = Ep
    save_Energies(j_save_E,2) = Ec
    save_Energies(j_save_E,3) = Ek
    	
end subroutine

subroutine Measure_of_Temperature()
implicit none
integer :: i
    v_rf_avg(:,1) = v_rf_avg(:,1) + v_x(:);
    v_rf_avg(:,2) = v_rf_avg(:,2) + v_y(:);
    v_rf_avg(:,3) = v_rf_avg(:,3) + v_z(:);
    S_rf_avg(:,1) = S_rf_avg(:,1) + S_rf_x
    S_rf_avg(:,2) = S_rf_avg(:,2) + S_rf_y

    do i = 1,3
        iRF(i) = iRF(i) + 1;
        if (iRF(i) == n_dt) then
		! print*, 'iRF =', iRF(1),iRF(2),iRF(3)
            iRF(i)           = 0;

            T_CM__LC(i)      = sum(v_rf_avg(:n_ions ,i) ) **2
            T_aux_LC(i)      = sum(v_rf_avg(:n_ions ,i)**2)
            S_aux_LC(i)		 = sum(S_rf_avg(:n_ions,i))
			
            v_rf_avg(:,i) = 0.0d0	
            S_rf_avg(:,i) = 0.0d0
			
        endif
    enddo
        
end subroutine

subroutine Check_save_temp()
implicit none
    if (n_save_temp == i_save_temp) then
        i_save_temp = 1;
        j_save_temp = j_save_temp + 1

        save_temperature(j_save_temp, 0 )   = dt*jj
		! print*, 'jj = ',jj, 'dt =', dt
		! print*, 'jj*dt =', jj*dt
        save_temperature(j_save_temp, 1 :3 ) = T_CM__LC
        save_temperature(j_save_temp, 4 :6 ) = T_aux_LC
        save_temperature(j_save_temp, 7 :8 ) = S_aux_LC(1:2)
        
        !dir$ if ( save_trj_cloud .eq. 1 )
		save_trj_Cloud(:,1,j_save_temp) =  sngl(r_x(:n_ions))
		save_trj_Cloud(:,2,j_save_temp) =  sngl(r_y(:n_ions))
		save_trj_Cloud(:,3,j_save_temp) =  sngl(r_z(:n_ions))
		save_trj_Cloud(:,4,j_save_temp) =  sngl(v_x(:n_ions))
		save_trj_Cloud(:,5,j_save_temp) =  sngl(v_y(:n_ions))
		save_trj_Cloud(:,6,j_save_temp) =  sngl(v_z(:n_ions))
		!dec$ endif
    else
        i_save_temp = i_save_temp +  1;
    endif
endsubroutine

subroutine Check_save_temp2()
implicit none
    if (n_save_temp == i_save_temp) then
        i_save_temp = 1;
        j_save_temp = j_save_temp + 1

        save_temperature(j_save_temp, 0 )   = dt*jj
		! print*, 'jj = ',jj, 'dt =', dt
		! print*, 'jj*dt =', jj*dt
        save_temperature(j_save_temp, 1 :3 ) = T_CM__LC
        save_temperature(j_save_temp, 4 :6 ) = T_aux_LC
        save_temperature(j_save_temp, 7 :8 ) = S_aux_LC(1:2)
        
        !dir$ if ( save_trj_cloud .eq. 1 )
		save_trj_Cloud(:,1,j_save_temp) =  sngl(r_x(:n_ions))
		save_trj_Cloud(:,2,j_save_temp) =  sngl(r_y(:n_ions))
		save_trj_Cloud(:,3,j_save_temp) =  sngl(r_z(:n_ions))
		save_trj_Cloud(:,4,j_save_temp) =  sngl(v_x(:n_ions))
		save_trj_Cloud(:,5,j_save_temp) =  sngl(v_y(:n_ions))
		save_trj_Cloud(:,6,j_save_temp) =  sngl(v_z(:n_ions))
		!dec$ endif
!~ 		print*, sum(T_aux_LC)*m_kb_x_inv_n_ions, T_threshold
		
        if (sum(T_aux_LC) > T_threshold) then
            flag_exit = 1;
        endif
    else
        i_save_temp = i_save_temp +  1;
    endif
endsubroutine

subroutine save_data()
implicit none
double precision, allocatable, dimension(:,:)   ::  save_final_position
integer :: i

    print*, trim(adjustl(str_file_xva))

    allocate(save_final_position(12,n_ions))

    save_final_position(1,:) = r_x ; save_final_position(2,:) = r_y ; save_final_position(3,:) = r_z
    save_final_position(4,:) = v_x ; save_final_position(5,:) = v_y ; save_final_position(6,:) = v_z
    save_final_position(7,:) = a1_x; save_final_position(8,:) = a1_y; save_final_position(9,:) = a1_z

    save_final_position(10,:)  = v_rf_avg(:,1);
    save_final_position(11,:)  = v_rf_avg(:,2);
    save_final_position(12,:)  = v_rf_avg(:,3);

    open(unit = 100, status='replace',file=trim(adjustl(str_file_xva))//'.bin',form='unformatted')  ! create a new file, or overwrite an existing one
        write(100) iRF
        write(100) save_final_position
    close(100)
    deallocate(save_final_position)

    open(101, file=trim(adjustl(str_file_xva))//'.info', status='replace', access='sequential', action='write')
        write(101,'(i16  , 3x, "%Number of LC ions")')                   n_ions
        write(101,'(i16  , 3x, "%Last index in the integration loop")')  j_start
        write(101,'(e16.9, 3x, "%Last time")')                           jj*dt
        write(101,'(i16  , 3x, "%Mass   of the LC ions")')               int(mass/amu)
        write(101,'(i16  , 3x, "%Charge of the LC ions")')               int(charge/qe)

        write(101,'(e16.9  , 3x, "%V_st[V]")')                             V_st
        write(101,'(e16.9  , 3x, "%V_rf[V]")')                             V_rf
        write(101,'(e16.9  , 3x, "%Omega/2pi[Hz]")')                       Omega/pi2
        write(101,'(e16.9  , 3x, "%r_0[m]")')                              r_0
        write(101,'(e16.9  , 3x, "%wz_LC/2 pi[Hz]")')                      wz_LC / pi2

        write(101,'(i16  , 3x, "%i_fly_Verlet")')                      i_fly_Verlet

    close(101)


!--------- Saving the Temperature Evolution -----------------------
!~     print*, 'Saving the Temperature Evolution in:',str_file_Temp_init
    open(9, file = trim(adjustl(str_file_Temp))//'.dat', status='replace', access='sequential', action='write')
    print*, 'j_save_temp Temp',j_save_temp
        do i = 1, j_save_temp
            write(9,221) save_temperature(i,0),&
                          m_kb_x_inv_n_ions2*save_temperature(i,1:3), &
                          m_kb_x_inv_n_ions *save_temperature(i,4:6), &
                          save_temperature(j_save_temp, 7 :8 )
        enddo

221     format( 9(1X,e27.19e3))
    close(9)

    deallocate(save_temperature)

!dir$ if(save_trj_cloud .eq. 1)
!--------- Saving the Cloud Evolution -----------------------
    open(unit = 111, status='replace',file=trim(adjustl(str_file_xva))//'_TrjN.bin',form='unformatted')  ! create a new file, or overwrite an existing one
    write(111) j_save_temp
    write(111) n_ions
    write(111) save_trj_Cloud
    close(111)
    deallocate(save_trj_Cloud)
!dir$ endif

endsubroutine

subroutine save_Energies_to_file()
implicit none
    ! Multiply by the correct cte before saving:
    save_Energies(:,2) = 0.5d+00*ke*charge**2 * save_Energies(:,2) ! Ec
    save_Energies(:,3) = (0.5d+00*mass)       * save_Energies(:,3) ! Ek
    print*, 'x',j_save_E
    open(unit = 11, status='replace',file='Energies_ndt100_vdt2.bin',form='unformatted')
        write(11) j_save_E
        write(11) size_Energy_array
        write(11) save_Energies
    close(11)
    print*, 'x',j_save_E
end subroutine save_Energies_to_file

subroutine manual_init_1ion()
implicit none
    iRF = 0
    j_start = 0 + 1 ! jj is the last one of the previous code, so the next one is jj+1
    print*,  iRF, j_start

    r_x  = 1
    r_y  = 0
    r_z  = 0

    v_x  = 0.0
    v_y  = 0.0
    v_z  = 0.0

    a1_x  = 0.0
    a1_y  = 0.0
    a1_z  = 0.0

    v_rf_avg(:,1) = 0.0
    v_rf_avg(:,2) = 0.0
    v_rf_avg(:,3) = 0.0
end 

subroutine init_from_file_Langevin()
implicit none
integer :: n_ions_aux
double precision, allocatable, dimension(:,:) :: xyz_temp
double precision :: dt_aux
print*, 'file to load', str_file_to_load
    open(unit = 10, status='old', file=trim(adjustl(str_file_to_load))//'.bin', form='unformatted')
        read(10) n_ions_aux
        read(10) jj
        read(10) dt_aux
        read(10) iRF
        allocate(xyz_temp(12,n_ions_aux))

        read(10) xyz_temp
    close(10)
    t_act   = jj*dt_aux
    t_next_n_dt = t_act + dt
    j_start = jj + 1 ! jj is the last one of the previous code, so the next one is jj+1
    print*, ' iRF, iRF, iRF, jstart '
    print*,  iRF, j_start

    r_x  = xyz_temp(1,:)
    r_y  = xyz_temp(2,:)
    r_z  = xyz_temp(3,:)

    v_x  = xyz_temp(4,:)
    v_y  = xyz_temp(5,:)
    v_z  = xyz_temp(6,:)

    a1_x  = xyz_temp(7,:)
    a1_y  = xyz_temp(8,:)
    a1_z  = xyz_temp(9,:)

    v_rf_avg(:,1) = xyz_temp(10,:)
    v_rf_avg(:,2) = xyz_temp(11,:)
    v_rf_avg(:,3) = xyz_temp(12,:)

    deallocate(xyz_temp)

end subroutine

subroutine create_files_names( )
implicit none
character(len=10 )  :: str_N,   str_N_aux
character(len=10 )  :: str_Vrf, str_Vrf_aux
character(len=19 )  :: str_Udc, str_Udc_aux
character(len=130)  :: str_file_aux, str_file_to_open
character(len=20 )  :: str_trap_aux


! N_ions
    if (n_ions<10) then
        write(str_N,"(I5.5)") n_ions
    elseif(n_ions < 100) then
        write(str_N,"(I5.5)") n_ions
    elseif(n_ions < 1000) then
        write(str_N,"(I5.5)") n_ions
    elseif(n_ions < 10000) then
        write(str_N,"(I5.5)") n_ions
    elseif(n_ions < 100000) then
        write(str_N,"(I5.5)") n_ions
    endif
    str_N_aux = '_N'//trim(adjustl(str_N))

! Vrf
    write(str_Vrf,"(I4.4)") int(V_rf)
    str_Vrf_aux = '_Vrf'//trim(adjustl(str_Vrf))

! Udc
    write(str_Udc,"(d13.4)") Udc
    str_Udc_aux = '_Udc'//trim(adjustl(str_Udc))//'V'

!! Udc
!    write(str_Udc,"(f6.0)") 1.d-3*wz_LC/pi2      sqrt(wz2(1))
!    str_Udc_aux = '_wz'//trim(adjustl(str_Udc))//'2pikHz'

    str_file_to_open = 'Langevin_cooling.info'
    open(15, file=trim(adjustl(str_file_to_open)), status='old', access='sequential', action='read')
		read(15,*)
		read(15,*)
		read(15,*)
		read(15,*)
		read(15,*)
		read(15,*)
		read(15,*)
		read(15,*)
		read(15,*)
		read(15,*)
		read(15,*)
        read(15,'(A)')    str_file_to_load
    close(15)
    str_trap_aux = 'SimuType0'

    ! Names used by the "from scrath" subroutine
    str_file_aux = trim(adjustl(str_trap_aux)) &
                // trim(adjustl(str_N_aux))    &
                // trim(adjustl(str_Vrf_aux))  &
                // trim(adjustl(str_Udc_aux))  &
                //trim(adjustl(str_extra))

    str_file_Temp = 'Temp_'//trim(adjustl(str_file_aux))
    str_file_xva  =  'xva_'//trim(adjustl(str_file_aux))
    str_file_trj  =  'trj_'//trim(adjustl(str_file_aux))

    print*, str_file_Temp
end subroutine

subroutine distribute_ions( )
implicit none
integer    :: i
integer                        :: nii, im
integer            :: n_by_core
integer, dimension(ni)   :: ia_tmp , ib_tmp
character(len=  10)  :: ia_str, ib_str
character(len=100)  :: str_ia, str_ib

    if (ni.ne.omp_get_max_threads ( )) then
        print*, 'ni (number of threads) is not the same in the makefile and in the code!'
        print*, ni, omp_get_max_threads ( )
        print*, 'Please, correct it, recompile and re-run'
        print*, 'Stoping code'
        stop
    endif

    i = modulo(n_ions, ni)
    n_by_core = floor(n_ions / real(ni))

    if (n_ions < ni) then
        nii = n_ions
        ia_tmp(:) = 1; ib_tmp(:) = 0; ! necessary for the ni > n_cut
        do im = 1, nii
            ia_tmp(im) = im
            ib_tmp(im) = im
        enddo
    else
        do im = 1, ni
            if (i == 0) then
                ia_tmp(im) = (im-1)*n_by_core+1
                ib_tmp(im) = ia_tmp(im) + n_by_core-1
            else
                if (im <= i) then
                    ia_tmp(im) = (im-1)*(n_by_core+1)+1
                    ib_tmp(im) = ia_tmp(im) + n_by_core
                else
                    ia_tmp(im) = n_by_core*(im - 1) + i + 1
                    ib_tmp(im) = ia_tmp(im) + n_by_core - 1
                endif
            endif
        enddo
    endif
    if (sum(ia - ia_tmp).ne.0) then
        print*, 'The ia, ib are not correct for the N and/or ni'
        print*, 'Update with the following: '

        str_ia = 'integer          , parameter :: ia(ni)=['
        str_ib = 'integer          , parameter :: ib(ni)=['

        do im = 1, ni
            if    (ia_tmp(im) < 10) then
                write(ia_str,"(I1.1)") ia_tmp(im)
            elseif (ia_tmp(im) < 100) then
                write(ia_str,"(I2.2)") ia_tmp(im)
            elseif (ia_tmp(im) < 1000) then
                write(ia_str,"(I3.3)") ia_tmp(im)
            else
                write(ia_str,"(I4.4)") ia_tmp(im)
            endif
            if (ib_tmp(im) < 10) then
                write(ib_str,"(I1.1)") ib_tmp(im)
            elseif (ib_tmp(im) < 100) then
                write(ib_str,"(I2.2)") ib_tmp(im)
            elseif (ib_tmp(im) < 1000) then
                write(ib_str,"(I3.3)") ib_tmp(im)
            else
                write(ib_str,"(I4.4)") ib_tmp(im)
            endif
            if (im == 1) then
                str_ia = trim(adjustl(str_ia))//trim(adjustl(ia_str))
                str_ib = trim(adjustl(str_ib))//trim(adjustl(ib_str))
            else
                str_ia = trim(adjustl(str_ia))//','//trim(adjustl(ia_str))
                str_ib = trim(adjustl(str_ib))//','//trim(adjustl(ib_str))
            endif
        enddo
        str_ia = trim(adjustl(str_ia))//'];'
        str_ib = trim(adjustl(str_ib))//'];'
        print*, str_ia
        print*, str_ib
        stop
    endif

endsubroutine

subroutine initialize_cos_sin() !Done
implicit none
integer :: i
do i = 1, n_dt
    cosi(i) = dcos((i-1)*dt_n)
    sini(i) = dsin((i-1)*dt_n)
enddo

end subroutine

end program main
