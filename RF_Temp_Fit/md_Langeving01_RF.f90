!~ Starting from md_v006_2D
!~ from md_Langevin_v73_1.f90

! I am using the vGB82 implemetation of the Langevin integrator
! ref: Skeel and Iwaguirrem Molecular Physics, v100, n24, p3885-3891  (2002)


!~ fx 360618.0563259394 Hz
!~ fy 360618.0563259394 Hz
!~ fz 240252.7346216562 Hz


include 'mkl_vsl.f90'
program main
!~ use ifport, only : random, seed
USE MKL_VSL_TYPE
USE MKL_VSL
use omp_lib
implicit none

!dir$ define potential_type = 2  ! 0 = Harmonique ; 1 = Polyfit; 2=Trigamma
!dir$ define simu_type = 3
!dir$ define integration_methode = 1  ! 0=Velocyt-Verlet; 1=Langeving_methode[vGB82]; 2=Langeving_methode[Impulse]

!! The Velocity-Verlet is not working well in the sense that I have only implemented viscosity, so there
!! is no source of heating....

include 'variable_definition.f90'

!*******************************************************************************
! For: vdRngGaussian
!*******************************************************************************
TYPE (VSL_STREAM_STATE) :: stream
integer brng,method
integer(kind=4) errcode

	print*, ''
    print*, '***********************************************************'
    print*, '        ==>   Welcome in md_Langevin01_RF.f90   <==        '
    print*, '                From Zero and Save to file                 '
    print*, '***********************************************************'

    call initialize_simu()
    time00 = omp_get_wtime()

    call free_fly()
    time_a = omp_get_wtime()
    print*, 'Total time', time_a - time00
    print*, 'Avg Time per step', (time_a - time00)/i_fly

    call create_files_names()
    print*, 'final jj', jj
	call save_info()
    call save_xva_to_file()
    call save_Temp_to_file()
    
    call save_trjectory1()
    call save_trjectoryN()
    
contains

subroutine free_fly()
implicit none
logical :: flag_check_performance = .True.
double precision :: time_per_step, expected_total_time

    j_end = j_end - 2 ! I need to remove 2 because of the zero and final steps happen outside the loop. 
    print*, 'start of free fly'
    print*,'jj',jj, 'jend',j_end
    call zero_time_step()


    j_aux = j_end - i_ss*n_dt
    do jj = 1, j_end
!~         print*, jj
        call update_positions()
        call update_accelerations()
        call update_velocities()
        call Measure_of_Temperature()
        call Check_save_temp()
        if (jj > j_aux) then
            save_trj1(jj - j_aux,1) = r_r(1)
            save_trj1(jj - j_aux,2) = r_y(1)
            save_trj1(jj - j_aux,3) = r_z(1)
            save_trj1(jj - j_aux,4) = v_r(1)
            save_trj1(jj - j_aux,5) = v_y(1)
            save_trj1(jj - j_aux,6) = v_z(1)
            save_trj1(jj - j_aux,7) = a1_r(1)
            save_trj1(jj - j_aux,8) = a1_y(1)
            save_trj1(jj - j_aux,9) = a1_z(1)

            save_trjN(jj - j_aux,:,1) = r_r(:)
            save_trjN(jj - j_aux,:,2) = r_y(:)
            save_trjN(jj - j_aux,:,3) = r_z(:)
            save_trjN(jj - j_aux,:,4) = v_r(:)
            save_trjN(jj - j_aux,:,5) = v_y(:)
            save_trjN(jj - j_aux,:,6) = v_z(:)
            save_trjN(jj - j_aux,:,7) = a1_r(:)
            save_trjN(jj - j_aux,:,8) = a1_y(:)
            save_trjN(jj - j_aux,:,9) = a1_z(:)
        endif
        if (flag_check_performance) then
            if ((jj-j_start+1) == 1000) then
                time_per_step = (omp_get_wtime() - time00) / 1000.0d0
                print*, 'Time for 1 step = ', time_per_step
                
                expected_total_time = (j_start + j_end - 2) * time_per_step
                print*, 'Expected Simulation Time[s]', expected_total_time
                print*, 'Expected Simulation Time[m]', expected_total_time/60.0
                print*, 'Expected Simulation Time[h]', expected_total_time/60.0/60.0
                flag_check_performance = .False.
            endif
        endif
    enddo
    call final_time_step()

    print*, 'jj after main loop and final step', jj
    save_trj1(jj - j_aux,1) = r_r(1)
    save_trj1(jj - j_aux,2) = r_y(1)
    save_trj1(jj - j_aux,3) = r_z(1)
    save_trj1(jj - j_aux,4) = v_r(1)
    save_trj1(jj - j_aux,5) = v_y(1)
    save_trj1(jj - j_aux,6) = v_z(1)
    save_trj1(jj - j_aux,7) = a1_r(1)
    save_trj1(jj - j_aux,8) = a1_y(1)
    save_trj1(jj - j_aux,9) = a1_z(1)

    save_trjN(jj - j_aux,:,1) = r_r(:)
    save_trjN(jj - j_aux,:,2) = r_y(:)
    save_trjN(jj - j_aux,:,3) = r_z(:)
    save_trjN(jj - j_aux,:,4) = v_r(:)
    save_trjN(jj - j_aux,:,5) = v_y(:)
    save_trjN(jj - j_aux,:,6) = v_z(:)
    save_trjN(jj - j_aux,:,7) = a1_r(:)
    save_trjN(jj - j_aux,:,8) = a1_y(:)
    save_trjN(jj - j_aux,:,9) = a1_z(:)


    print*, 'i_ss*n_dt+10',i_ss*n_dt+10
    print*, 'jj - j_aux',jj - j_aux
!~     jj = jj + 1

endsubroutine

subroutine save_trjectory1()
implicit none
integer :: i
    open(unit = 11, status='replace',file='trj_L.bin',form='unformatted')
        write(11) jj-j_aux
        write(11) save_trj1(:jj-j_aux,:)
    close(11)
end subroutine save_trjectory1

subroutine save_trjectoryN()
implicit none
integer :: i
!~     print*, save_trjN(jj-j_aux,:,e:)
    print*,'hi'
    open(unit = 11, status='replace',file='trjN_L.bin',form='unformatted')
        write(11) n_ions
        write(11) jj-j_aux
        write(11) i_ss*n_dt+1
        write(11) save_trjN
    close(11)
end subroutine save_trjectoryN

subroutine zero_time_step()
implicit none

    alphaL  = dsqrt(cte1)

    errcode=vdrnggaussian( method, stream, n_ions,Z0_r, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z0_z, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z0_y, 0.0d0, 1.0d0)

    v_z  = e_eta_mass2 * (v_z + dt*w_p*a2_z + alphaL*Z0_z) ! v_{1/2}
    v_r  = e_eta_mass2 * (v_r + dt*w_p*a2_r + alphaL*Z0_r) ! v_{1/2}
    v_y  = e_eta_mass2 * (v_y + dt*w_p*a2_y + alphaL*Z0_y) ! v_{1/2}
end subroutine

subroutine update_positions()
implicit none
integer :: im
!-----------------------------------------------------------------------
!   One time step
    !$omp parallel default(none) private(im) shared (r_z,v_z,r_r,v_r,r_y,v_y) firstprivate(ia,ib)
    im = omp_get_thread_num()+1
    ! This falls into the "fused multiply add" category: U = U + B*X
    ! It could be speed up!
    r_r(ia(im): ib(im)) = r_r(ia(im): ib(im)) + v_r(ia(im): ib(im))*e_eta_mass3;
    r_z(ia(im): ib(im)) = r_z(ia(im): ib(im)) + v_z(ia(im): ib(im))*e_eta_mass3;
    r_y(ia(im): ib(im)) = r_y(ia(im): ib(im)) + v_y(ia(im): ib(im))*e_eta_mass3;
    !$omp end parallel
end subroutine

subroutine update_accelerations()
implicit none
!---------------  Locals  ----------------------------------------------
integer :: i, im,j
double precision :: rji(3)
double precision :: r2inv
double precision, parameter :: softening = 1.00D-20
double precision    :: cte_aux(2)
double precision    :: cos_jj

    cos_jj = cosi(mod(jj,n_dt)+1)    
    cte_aux(1) = a_cte_Quad_RF_LC(1) + a_cte_Quad_RF_LC(3)*cos_jj
    cte_aux(2) = a_cte_Quad_RF_LC(2) - a_cte_Quad_RF_LC(3)*cos_jj
! -------------------------------------------------
!~     print*, 'RF Langevin', a_cte_Quad_RF_LC
!~     pause
    !$omp parallel default(none) &
    !$omp private(im, i,j,rji,r2inv) &
    !$omp firstprivate(ia, ib, r_z, r_r,r_y,cte_aux) &
    !$omp shared(a2_z,a2_r,a2_y)
    im = omp_get_thread_num() + 1
! RF contribution to acceleration
    a2_r(ia(im):ib(im)) =          cte_aux(1) * r_r(ia(im):ib(im))
    a2_y(ia(im):ib(im)) =          cte_aux(2) * r_y(ia(im):ib(im))
    a2_z(ia(im):ib(im)) = a_cte_Quad_RF_LC(4) * r_z(ia(im):ib(im))
! Pseudo-pot contribution to acceleration
    !    a2_r(ia(im): ib(im)) = - wr2*r_r(ia(im): ib(im))
    !    a2_y(ia(im): ib(im)) = - wy2*r_y(ia(im): ib(im))
    !    a2_z(ia(im): ib(im)) = - wz2*r_z(ia(im): ib(im))
! Coulombian repulsion see a_Coulomb_1sp_2sp()
    do i = ia(im), ib(im)
        do j = 1, n_ions
            rji(1)  = r_z(j) - r_z(i)
            rji(2)  = r_r(j) - r_r(i)
            rji(3)  = r_y(j) - r_y(i)
            r2inv   = 1.00D+00/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3) + softening)
            r2inv   = r2inv * r2inv * r2inv * alpha
            a2_z(i) = a2_z(i) - rji(1)*r2inv
            a2_r(i) = a2_r(i) - rji(2)*r2inv
            a2_y(i) = a2_y(i) - rji(3)*r2inv
        enddo
    enddo
    !$omp end parallel
!~     print*, cos_jj
!~     print*,'a'
!~     print*, a2_r(1), a2_y(1), a2_z(1)
!~     print*, r_r(1), r_y(1), r_z(1)
!~     print*, v_r(1), v_y(1), v_z(1)
!~         pause

endsubroutine

subroutine update_velocities()
implicit none
integer :: im
real ( kind = 8 ), parameter :: A1 = e_eta_mass2*e_eta_mass2
real ( kind = 8 ), parameter :: A2 = e_eta_mass2*dt_S2
real ( kind = 8 ), parameter :: A3 = -e_eta_mass2*dt_S
real ( kind = 8 ) :: A4, A5
    betaL   = cte2/alphaL             ; A4 = e_eta_mass2*betaL
    alphaL  = dsqrt(cte13-betaL*betaL); A5 = e_eta_mass2*alphaL
    
    errcode=vdrnggaussian( method, stream, n_ions,Z1_r, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z1_z, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z1_y, 0.0d0, 1.0d0)

    !$omp parallel default(none) private(im) shared (v_z,a1_z,a2_z,v_r,a1_r,a2_r,v_y,a1_y,a2_y, Z1_r, Z0_r,Z1_y, Z0_y, Z1_z, Z0_z) firstprivate(ia,ib,A4,A5)
    im = omp_get_thread_num()+1

    v_r(ia(im):ib(im))  = A1*v_r(ia(im):ib(im)) + A2*a2_r(ia(im):ib(im)) + A3*a1_r(ia(im):ib(im)) + A4*Z0_r(ia(im):ib(im)) + A5*Z1_r(ia(im):ib(im));
    a1_r(ia(im):ib(im)) = a2_r(ia(im):ib(im)); 
    Z0_r(ia(im):ib(im)) = Z1_r(ia(im):ib(im));

    v_z(ia(im):ib(im))  = A1*v_z(ia(im):ib(im)) + A2*a2_z(ia(im):ib(im)) + A3*a1_z(ia(im):ib(im)) + A4*Z0_z(ia(im):ib(im)) + A5*Z1_z(ia(im):ib(im));
    a1_z(ia(im):ib(im)) = a2_z(ia(im):ib(im)); 
    Z0_z(ia(im):ib(im)) = Z1_z(ia(im):ib(im));

    v_y(ia(im):ib(im))  = A1*v_y(ia(im):ib(im)) + A2*a2_y(ia(im):ib(im)) + A3*a1_y(ia(im):ib(im)) + A4*Z0_y(ia(im):ib(im)) + A5*Z1_y(ia(im):ib(im));
    a1_y(ia(im):ib(im)) = a2_y(ia(im):ib(im)); 
    Z0_y(ia(im):ib(im)) = Z1_y(ia(im):ib(im));
    !$omp end parallel
end subroutine

subroutine Measure_of_Temperature()
implicit none
integer :: i
    v_rf_avg(:,1) = v_rf_avg(:,1) + v_r(:);
    v_rf_avg(:,2) = v_rf_avg(:,2) + v_y(:);
    v_rf_avg(:,3) = v_rf_avg(:,3) + v_z(:);

    do i = 1,3
        iRF(i) = iRF(i) + 1;
        if (iRF(i) == n_dt) then
!~ 		print*, 'iRF =', iRF(1),iRF(2),iRF(3)
            iRF(i)        = 0;
            T_CM (i)      = sum(v_rf_avg(:, i) )**2
            T_aux(i)      = sum(v_rf_avg(:, i  )**2)
            v_rf_avg(:,i) = 0.0d0
        endif
    enddo
    
    return
end subroutine

subroutine Check_save_temp()
implicit none
    if (n_save_temp == i_save_temp) then
        i_save_temp = 1;
        j_save_temp = j_save_temp + 1

        save_temperature(j_save_temp, 0    ) = jj*dt
		! print*, 'jj = ',jj, 'dt =', dt
		! print*, 'jj*dt =', jj*dt
        save_temperature(j_save_temp, 1 :3 ) = T_CM
        save_temperature(j_save_temp, 4 :6 ) = T_aux
    else
        i_save_temp = i_save_temp + 1;
    endif
endsubroutine

subroutine final_time_step()
implicit none
    call update_positions();
    call update_accelerations()
    call update_velocities_final_step()
end subroutine

subroutine update_velocities_final_step()
implicit none
integer :: im
!~ beta    = b/alpha
!~ alpha   = sqrt(c-beta**2)
!~ Z1      = random.normal()
!~ v_0_12  = e_gamma2 * v_0_12 + dt*w_m*F_1 + dt*S*(F_1-F_0) + beta*Z0 + alpha*Z1
    betaL   = cte2/alphaL
!~     print*, cte3-betaL*betaL
    alphaL  = dsqrt(abs(cte3-betaL*betaL))
    errcode=vdrnggaussian( method, stream, n_ions,Z1_r, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z1_z, 0.0d0, 1.0d0)
    errcode=vdrnggaussian( method, stream, n_ions,Z1_y, 0.0d0, 1.0d0)

    !$omp parallel default(none) private(im) shared (v_z,a1_z,a2_z,v_r,a1_r,a2_r, v_y,a1_y,a2_y, Z1_r, Z0_r, Z1_y, Z0_y,Z1_z, Z0_z) firstprivate(ia,ib, betaL, alphaL)
    im = omp_get_thread_num()+1
!dir$ if (integration_methode.eq.0) ! Velocity-Verlet

!dir$ elseif (integration_methode.eq.1) ! Langeving_methode
    v_r(ia(im):ib(im)) = e_eta_mass2*v_r(ia(im):ib(im)) + dt*w_m*a2_r(ia(im):ib(im)) + dt_S*(a2_r(ia(im):ib(im))-a1_r(ia(im):ib(im))) + betaL*Z0_r(ia(im):ib(im)) + alphaL*Z1_r(ia(im):ib(im));
    v_z(ia(im):ib(im)) = e_eta_mass2*v_z(ia(im):ib(im)) + dt*w_m*a2_z(ia(im):ib(im)) + dt_S*(a2_z(ia(im):ib(im))-a1_z(ia(im):ib(im))) + betaL*Z0_z(ia(im):ib(im)) + alphaL*Z1_z(ia(im):ib(im));
    v_y(ia(im):ib(im)) = e_eta_mass2*v_y(ia(im):ib(im)) + dt*w_m*a2_y(ia(im):ib(im)) + dt_S*(a2_y(ia(im):ib(im))-a1_y(ia(im):ib(im))) + betaL*Z0_y(ia(im):ib(im)) + alphaL*Z1_y(ia(im):ib(im));
!dir$ elseif (integration_methode.eq.2) ! langeving_methode: Impulse
    v_r(ia(im):ib(im)) = e_eta_mass2*v_r(ia(im):ib(im)) + dt*w_m*a2_r(ia(im):ib(im)) + betaL*Z0_r(ia(im):ib(im)) + alphaL*Z1_r(ia(im):ib(im));
    v_z(ia(im):ib(im)) = e_eta_mass2*v_z(ia(im):ib(im)) + dt*w_m*a2_z(ia(im):ib(im)) + betaL*Z0_z(ia(im):ib(im)) + alphaL*Z1_z(ia(im):ib(im));
    v_y(ia(im):ib(im)) = e_eta_mass2*v_y(ia(im):ib(im)) + dt*w_m*a2_y(ia(im):ib(im)) + betaL*Z0_y(ia(im):ib(im)) + alphaL*Z1_y(ia(im):ib(im));
!dir$ endif
    !$omp end parallel
end subroutine

subroutine initialize_simu() !Done
implicit none
    
    ! First of all: check that cte1 > 0
    if ( cte1 < 0 ) then
        print*, 'cte1 is negative!'
        print*, 'cte1 = ', cte1
        print*, 'Stopping code.'
        stop
    endif
    
    call init_random_generator()
    call distribute_ions()
    call initialize_ions()

    call initialize_cos()

    
    iRF(1) = 0; 
    iRF(2) = floor(n_dt/2.0); 
    iRF(3) = 0;
    v_rf_avg = 0.0d0

    j_end   = i_fly*n_dt - modulo(i_fly*n_dt,n_dt) + int(0.25*n_dt)
    print*,  'j_end',j_end

    jj = int(j_end/real(n_save_temp))
    allocate(save_temperature(jj,0:6))
    
    allocate(save_trj1(i_ss*n_dt+10,9))
    allocate(save_trjN(i_ss*n_dt+1,n_ions,9))
    
    call update_accelerations()

    a1_r = a2_r; a2_r = 0.0d0;
    a1_y = a2_y; a2_y = 0.0d0;
    a1_z = a2_z; a2_z = 0.0d0;
    

end subroutine

subroutine initialize_cos() !Done
implicit none
integer :: i
do i = 1, n_dt
    cosi(i) = dcos((i-1)*dt_n)
enddo

end subroutine

subroutine initialize_ions() !Done
implicit none
double precision, parameter :: l0(3) = [10.0e-6,10.0e-6,10.0e-6];

    errcode=vdrnggaussian( method, stream, n_ions,r_r, 0.0d0, l0(1))
    errcode=vdrnggaussian( method, stream, n_ions,r_y, 0.0d0, l0(3))
    errcode=vdrnggaussian( method, stream, n_ions,r_z, 0.0d0, l0(2))
    
!~     r_r(1) = 1e-6
!~     r_y(1) = 2e-6
!~     r_z(1) = 3e-6

    v_z  = 0.0d0; a1_z = 0.0d0; a2_z = 0.0d0;
    v_r  = 0.0d0; a1_r = 0.0d0; a2_r = 0.0d0;
    v_y  = 0.0d0; a1_y = 0.0d0; a2_y = 0.0d0;

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

subroutine create_files_names( )
implicit none
character(len=10 )  :: str_N  ,   str_N_aux
character(len=130)  :: str_trap_aux
character(len=50 )  :: str_T,str_F

    write(str_N,"(I4.4)") n_ions
    str_N_aux = '_N'//trim(adjustl(str_N))

    str_trap_aux = '3D_Harmo'


    ! Temperature
    if     (Temperature .lt. 1.0E-06) then
        write(str_T,"(I3.3)") int(Temperature*1.0d9)
        str_T = '_T'//trim(adjustl(str_T))//'nK'
    elseif (Temperature .lt. 1.0E-03) then
        write(str_T,"(I3.3)") int(Temperature*1.0d6)
        str_T = '_T'//trim(adjustl(str_T))//'uK'
    elseif (Temperature .lt. 1.0E-00) then
        write(str_T,"(I3.3)") int(Temperature*1.0d3)
        str_T = '_T'//trim(adjustl(str_T))//'mK'
    endif
    ! Friction coefficient
    write(str_F,"(d9.2)") abs(eta)
    str_F = '_F'//trim(adjustl(str_F))//'Kg_s'

    ! Names used by the "from scrath" subroutine
    str_file_aux = trim(adjustl(str_trap_aux)) &
                // trim(adjustl(str_N_aux))    &
                // trim(adjustl(str_T))    &
                // trim(adjustl(str_F))    &
                // trim(adjustl(str_extra_Lan))

!~ print*, str_file_aux
!~ stop
end subroutine

subroutine save_Temp_to_file()
implicit none
character(len=130)  :: str_file_trj_tmp
double precision, parameter :: m_kb_x_inv_n_ions  = mass/(kb*real(n_ions)*n_dt**2   )
double precision, parameter :: m_kb_x_inv_n_ions2 = mass/(kb*real(n_ions)**2*n_dt**2)

    str_file_trj_tmp = 'Temp_'//trim(adjustl(str_file_aux))//'.bin'
    open(unit = 11, status='replace',file=trim(str_file_trj_tmp),form='unformatted')  ! create a new file, or overwrite an existing one
        write(11) n_ions
        write(11) j_save_temp
        write(11) dt*n_save_temp
        write(11) eta
        write(11) Temperature

        write(11) save_temperature(:j_save_temp,0),&
                  save_temperature(:j_save_temp,1:3)*m_kb_x_inv_n_ions2,&
                  save_temperature(:j_save_temp,4:6)*m_kb_x_inv_n_ions
                  
    close(11)

!~     print*, save_temperature(1:10,1)*m_kb_x_inv_n_ions2
!~     print*, save_temperature(1:10,2)*m_kb_x_inv_n_ions2


    deallocate(save_temperature)
end subroutine

subroutine save_xva_to_file()
implicit none
character(len=130)  :: str_file_trj_tmp
double precision, allocatable, dimension(:,:)   ::  save_final_position

    str_file_trj_tmp = 'xva_'//trim(adjustl(str_file_aux))//'.bin'
    print*, 'Saving Final positions to:', trim(str_file_trj_tmp)
   
    allocate(save_final_position(12,n_ions))

    save_final_position(1 ,:) = r_r ; save_final_position(2,:) = r_y ; save_final_position(3,:) = r_z
    save_final_position(4 ,:) = v_r ; save_final_position(5,:) = v_y ; save_final_position(6,:) = v_z
    save_final_position(7 ,:) = a1_r; save_final_position(8,:) = a1_y; save_final_position(9,:) = a1_z

    save_final_position(10,:) = v_rf_avg(:,1);
    save_final_position(11,:) = v_rf_avg(:,2);
    save_final_position(12,:) = v_rf_avg(:,3);

   
    open(unit = 11, status='replace',file=trim(str_file_trj_tmp),form='unformatted')  ! create a new file, or overwrite an existing one
        write(11) n_ions
        write(11) jj
        write(11) dt
        write(11) iRF
        write(11) save_final_position
    close(11)
    deallocate(save_final_position)
end subroutine

subroutine save_info()
implicit none
character(len=130)  :: str_file_trj_tmp, str_file_info_tmp
double precision, allocatable, dimension(:,:)   ::  save_final_position

	str_file_trj_tmp = 'xva_'//trim(adjustl(str_file_aux))
    str_file_info_tmp = 'Langevin_cooling.info'

	open(unit = 10, file=trim(adjustl(str_file_info_tmp)), status='replace', access='sequential', action='write')
        write(10,'(i16  , 3x, "%Number of LC ions")')                   n_ions
        write(10,'(i16  , 3x, "%Last index in the integration loop")')  jj
        write(10,'(e16.9, 3x, "%Last time")')                           jj*dt
        write(10,'(i16  , 3x, "%Mass   of the LC ions")')               int(mass/amu)
        write(10,'(i16  , 3x, "%Charge of the LC ions")')               int(charge/qe)

        write(10,'(e16.9  , 3x, "%V_st[V]")')                             V_st
		write(10,'(e16.9  , 3x, "%U_dc[V]")')                             Udc
        write(10,'(e16.9  , 3x, "%V_rf[V]")')                             V_rf
        write(10,'(e16.9  , 3x, "%Omega/2pi[Hz]")')                       Omega/pi2

		write(10,'(e16.9  , 3x, "%Eta")')                            Eta
        write(10,'(i16  , 3x, "%i_free__fly")')                      i_fly
        write(10,'(A)') trim(adjustl(str_file_trj_tmp))
!~         write(10,'(i16  , 3x, "%i_laser_fly")')                      '##v##ADRIEN##v##'
!~         write(10,'(e16.9  , 3x, "%detuning")')                       detuning
!~         write(10,'(e16.9  , 3x, "%saturation")')                     saturation

    close(10)
end subroutine save_info

subroutine init_random_generator()
implicit none
    brng=VSL_BRNG_MCG31
    method=VSL_RNG_METHOD_GAUSSIAN_ICDF

    CALL RANDOM_SEED()
    call random_number(rand_num)
    rand_seed = int(rand_num(1)*1e4)
    errcode=vslnewstream( stream, brng, rand_seed )
    
    return
end subroutine init_random_generator


end program main
