! -*- f90 -*-
subroutine solveSAS(J_fullstep, Q_fullstep, SAS_args, P_list, weights_fullstep, sT_init_fullstep, dt, &
                    verbose, debug, warning, jacobian, &
                    mT_init_fullstep, C_J_fullstep, alpha_fullstep, k1_fullstep, C_eq_fullstep, C_old, &
                    n_substeps, component_type, numcomponent_list, numargs_list, numflux, numsol, max_age, &
                    timeseries_length, output_these_fullsteps, num_output_fullsteps, numcomponent_total, numargs_total, &
                    num_scheme, sT_outputstep, pQ_outputstep, WaterBalance_outputstep, &
                    mT_outputstep, mQ_outputstep, mR_outputstep, C_Q_fullstep, ds_outputstep, dm_outputstep, &
                    dC_fullstep, SoluteBalance_outputstep)
   use cdf_gamma_mod
   use cdf_beta_mod
   !use cdf_normal_mod
   implicit none

   ! Start by declaring and initializing all the variables we will be using
   integer, intent(in) :: n_substeps, numflux, numsol, max_age, num_scheme, &
                          timeseries_length, numcomponent_total, numargs_total, num_output_fullsteps
   real(8), intent(in) :: dt
   logical, intent(in) :: verbose, debug, warning, jacobian
   real(8), intent(in), dimension(0:timeseries_length - 1) :: J_fullstep
   real(8), intent(in), dimension(0:timeseries_length - 1, 0:numflux - 1) :: Q_fullstep
   real(8), intent(in), dimension(0:timeseries_length - 1, 0:numcomponent_total - 1) :: weights_fullstep
   real(8), intent(in), dimension(0:numargs_total - 1, 0:timeseries_length - 1) :: SAS_args
   real(8), intent(in), dimension(0:numargs_total - 1, 0:timeseries_length - 1) :: P_list
   real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: C_J_fullstep
   real(8), intent(in), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: alpha_fullstep
   real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: k1_fullstep
   real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: C_eq_fullstep
   real(8), intent(in), dimension(0:numsol - 1) :: C_old
   real(8), intent(in), dimension(0:max_age - 1) :: sT_init_fullstep
   real(8), intent(in), dimension(0:max_age - 1, 0:numsol - 1) :: mT_init_fullstep
   integer, intent(in), dimension(0:num_output_fullsteps - 1) :: output_these_fullsteps
   integer, intent(in), dimension(0:numcomponent_total - 1) :: component_type
   integer, intent(in), dimension(0:numflux - 1) :: numcomponent_list
   integer, intent(in), dimension(0:numcomponent_total - 1) :: numargs_list
   real(8), intent(out), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: C_Q_fullstep
   real(8), intent(out), dimension(0:timeseries_length - 1, 0:numargs_total - 1, 0:numflux - 1, 0:numsol - 1) :: dC_fullstep
   real(8), intent(out), dimension(0:num_output_fullsteps, 0:max_age - 1) :: sT_outputstep
   real(8), intent(out), dimension(0:num_output_fullsteps, 0:numsol - 1, 0:max_age - 1) :: mT_outputstep
   real(8), intent(out), dimension(0:num_output_fullsteps, 0:numargs_total - 1, 0:max_age - 1) :: ds_outputstep
   real(8), intent(out), dimension(0:num_output_fullsteps, 0:numargs_total - 1, 0:numsol - 1, 0:max_age - 1) :: dm_outputstep
   real(8), intent(out), dimension(0:num_output_fullsteps - 1, 0:numflux - 1, 0:max_age - 1) :: pQ_outputstep
   real(8), intent(out), dimension(0:num_output_fullsteps - 1, 0:numflux - 1, 0:numsol - 1, 0:max_age - 1) :: mQ_outputstep
   real(8), intent(out), dimension(0:num_output_fullsteps - 1, 0:numsol - 1, 0:max_age - 1) :: mR_outputstep
   real(8), intent(out), dimension(0:num_output_fullsteps - 1, 0:max_age - 1) :: WaterBalance_outputstep
   real(8), intent(out), dimension(0:num_output_fullsteps - 1, 0:numsol - 1, 0:max_age - 1) :: SoluteBalance_outputstep
   real(8), dimension(0:timeseries_length - 1, 0:numargs_total - 1, 0:numflux - 1) :: dW_outputstep
   real(8), dimension(0:timeseries_length - 1, 0:numflux - 1) :: P_old_fullstep
   integer, dimension(0:numcomponent_total) :: args_index_list
   integer, dimension(0:numflux) :: component_index_list
   real(8), dimension(0:timeseries_length*n_substeps, 0:1) :: STcum_topbot_start
   integer, dimension(0:timeseries_length*n_substeps - 1, 0:numcomponent_total - 1, 0:1) :: leftbreakpt_topbot
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1) :: pQ_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1) :: pQ_aver
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_aver
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mR_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mR_aver
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1) :: fs_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1) :: fs_aver
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numflux - 1) :: fsQ_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numflux - 1) :: fsQ_aver
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numsol - 1) :: fm_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numsol - 1) :: fm_aver
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numflux - 1, 0:numsol - 1) :: fmQ_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numflux - 1, 0:numsol - 1) :: fmQ_aver
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numsol - 1) :: fmR_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numsol - 1) :: fmR_aver
   real(8), dimension(0:timeseries_length*n_substeps - 1) :: sT_start
   real(8), dimension(0:timeseries_length*n_substeps - 1) :: sT_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mT_start
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mT_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1) :: ds_start
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1) :: ds_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numsol - 1) :: dm_start
   real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numargs_total - 1, 0:numsol - 1) :: dm_temp
   real(8), dimension(0:timeseries_length*n_substeps - 1) :: STcum_in
   integer, dimension(0:timeseries_length*n_substeps - 1) :: jt_fullstep_at_
   integer, dimension(0:timeseries_length*n_substeps - 1) :: jt_substep_at_
   real(8), dimension(0:numargs_total - 1, 0:timeseries_length - 1) :: grad_precalc
   integer :: substep, iT_fullstep, iT_substep, iT_prev_fullstep, jt_is_which_substep, jt_fullstep
   real(8) :: one8, norm
   real(8) :: dS, dP, dSe, dPe, dSs, dPs
   real(8) :: dt_substep, dt_numerical_solution
   real(8), dimension(2) :: rk2_coeff
   real(8), dimension(3) :: rk2_stepfraction
   real(8), dimension(4) :: rk4_coeff
   real(8), dimension(5) :: rk4_stepfraction
   character(len=128) :: tempdebugstring
   integer :: iq, s, total_num_substeps, ip, ic, c, rk, outputstep, jt_c
   integer :: carry
   integer :: leftbreakpt
   integer :: ia
   integer :: i
   real(8) :: start, finish

   C_Q_fullstep = 0.
   sT_outputstep = 0.
   mT_outputstep = 0.
   ds_outputstep = 0.
   dm_outputstep = 0.
   dC_fullstep = 0.
   dW_outputstep = 0.
   pQ_outputstep = 0.
   mQ_outputstep = 0.
   mR_outputstep = 0.
   WaterBalance_outputstep = 0.
   SoluteBalance_outputstep = 0.
   P_old_fullstep = 1.
   args_index_list = 0
   component_index_list = 0
   STcum_topbot_start = 0.
   leftbreakpt_topbot = 0
   pQ_temp = 0.
   pQ_aver = 0.
   mQ_temp = 0.
   mQ_aver = 0.
   mR_temp = 0.
   mR_aver = 0.
   fs_temp = 0.
   fs_aver = 0.
   fsQ_temp = 0.
   fsQ_aver = 0.
   fm_temp = 0.
   fm_aver = 0.
   fmQ_temp = 0.
   fmQ_aver = 0.
   fmR_temp = 0.
   fmR_aver = 0.
   sT_start = 0.
   sT_temp = 0.
   mT_start = 0.
   mT_temp = 0.
   ds_start = 0.
   ds_temp = 0.
   dm_start = 0.
   dm_temp = 0.
   iT_prev_fullstep = -1

   call f_verbose('...Initializing arrays...')
   one8 = 1.0
   rk4_stepfraction = (/0.0D0, 0.5D0, 0.5D0, 1.0D0, 1.0D0/)
   rk4_coeff = (/1./6, 2./6, 2./6, 1./6/)
   rk2_stepfraction = (/0.0D0, 1.0D0, 1.0D0/)
   rk2_coeff = (/1./2, 1./2/)
   norm = 1.0/n_substeps/n_substeps

   ! The list of probabilities in each sas function is a 1-D array.
   ! args_index_list gives the starting index of the probabilities (P) associated
   ! with each flux
   args_index_list(0) = 0
   component_index_list(0) = 0
   do iq = 0, numflux - 1
      component_index_list(iq + 1) = component_index_list(iq) + numcomponent_list(iq)
      do ic = component_index_list(iq), component_index_list(iq + 1) - 1
         args_index_list(ic + 1) = args_index_list(ic) + numargs_list(ic)
      end do
   end do
   call f_debug('args_index_list', one8*args_index_list(:))

   do iq = 0, numflux - 1
      do ic = component_index_list(iq), component_index_list(iq + 1) - 1
         if (component_type(ic) == -1) then
            do ia = 0, numargs_list(ic) - 1
               grad_precalc(args_index_list(ic) + ia, :) = &
                  (P_list(args_index_list(ic) + ia + 1, :) - P_list(args_index_list(ic) + ia, :)) &
                  /(SAS_args(args_index_list(ic) + ia + 1, :) - SAS_args(args_index_list(ic) + ia, :))
            end do
         end if
      end do
   end do

   ! modify the number of ages and the timestep by a facotr of n_substeps
   total_num_substeps = timeseries_length*n_substeps
   dt_substep = dt/n_substeps

   call f_verbose('...Setting initial conditions...')
   sT_outputstep(0, :) = sT_init_fullstep
   mT_outputstep(0, :, :) = mT_init_fullstep

   call f_verbose('...Starting main loop...')
   ! Loop over ages
   do iT_fullstep = 0, max_age - 1

      ! Start the substep loop
      do substep = 0, n_substeps - 1

         iT_substep = iT_fullstep*n_substeps + substep

         ! jt_fullstep_at_(c) maps characteristic index c to the full timestep it is currently intersecting
         do c = 0, total_num_substeps - 1
            jt_substep_at_(c) = mod(c + iT_substep, total_num_substeps)
            jt_is_which_substep = mod(jt_substep_at_(c), n_substeps)
            jt_fullstep_at_(c) = (jt_substep_at_(c) - jt_is_which_substep)/n_substeps
         end do

         pQ_aver = 0
         mQ_aver = 0
         mR_aver = 0

         ! Apply the initial conditions
         if (iT_substep > 0) then
            sT_start(total_num_substeps - iT_substep) = sT_init_fullstep(iT_prev_fullstep)
            mT_start(total_num_substeps - iT_substep, :) = mT_init_fullstep(iT_prev_fullstep, :)
         end if

         sT_temp = sT_start
         mT_temp = mT_start

         call f_debug('sT_temp                           ', sT_temp)

         !if (jacobian) then
         !fs_aver = 0
         !fsQ_aver = 0
         !fm_aver = 0
         !fmQ_aver = 0
         !fmR_aver = 0
         !if (iT_substep > 0) then
         !ds_start(total_num_substeps - iT_substep, :) = 0.
         !dm_start(total_num_substeps - iT_substep, :, :) = 0.
         !end if
         !ds_temp = ds_start
         !dm_temp = dm_start
         !end if

         select case (num_scheme)
         case (1)
            ! This is the forward Euler
            call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, one8)
            if (iT_substep == 0) then
               call new_state(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, one8)
               call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, one8)
            end if
            call add_to_average(pQ_temp, mQ_temp, mR_temp, one8)
            call new_state(sT_temp, mT_temp, pQ_aver, mQ_aver, mR_aver, one8)

         case (2)
            ! This is the Runge-Kutta 2nd order algorithm
            rk = 1
            call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk2_stepfraction(rk))
            call add_to_average(pQ_temp, mQ_temp, mR_temp, rk2_coeff(rk))
            rk = 2
            call new_state(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk2_stepfraction(rk))
            call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk2_stepfraction(rk))
            call add_to_average(pQ_temp, mQ_temp, mR_temp, rk2_coeff(rk))
            rk = 3
            call new_state(sT_temp, mT_temp, pQ_aver, mQ_aver, mR_aver, rk2_stepfraction(rk))

         case (4)
            ! This is the Runge-Kutta 4th order algorithm
            rk = 1
            call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk4_stepfraction(rk))
            call add_to_average(pQ_temp, mQ_temp, mR_temp, rk4_coeff(rk))
            rk = 2
            call new_state(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk4_stepfraction(rk))
            call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk4_stepfraction(rk))
            call add_to_average(pQ_temp, mQ_temp, mR_temp, rk4_coeff(rk))
            rk = 3
            call new_state(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk4_stepfraction(rk))
            call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk4_stepfraction(rk))
            call add_to_average(pQ_temp, mQ_temp, mR_temp, rk4_coeff(rk))
            rk = 4
            call new_state(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk4_stepfraction(rk))
            call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk4_stepfraction(rk))
            call add_to_average(pQ_temp, mQ_temp, mR_temp, rk4_coeff(rk))
            rk = 5
            call new_state(sT_temp, mT_temp, pQ_aver, mQ_aver, mR_aver, rk4_stepfraction(rk))
         end select

         !! This is the Runge-Kutta 4th order algorithm
         !do rk = 1, 5
         !if (rk > 1) then
         !call new_state(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk_stepfraction(rk))
         !end if
         !if (rk < 5) then
         !call get_flux(sT_temp, mT_temp, pQ_temp, mQ_temp, mR_temp, rk_stepfraction(rk))
         !call add_to_average(pQ_temp, mQ_temp, mR_temp, rk_coeff(rk))
         !end if
         !if (rk == 4) then
               !!call f_debug('FINALIZE  rk           ', (/rk*one8, iT_substep*one8/))
         !pQ_temp = pQ_aver
         !mQ_temp = mQ_aver
         !mR_temp = mR_aver
               !!if (jacobian) then
               !!fs_temp = fs_aver
               !!fsQ_temp = fsQ_aver
               !!fm_temp = fm_aver
               !!fmR_temp = fmR_aver
               !!fmQ_temp = fmQ_aver
               !!end if
         !call f_debug('pQ_aver                ', pQ_aver(:, 0))
         !end if
         !end do
         !call f_debug_blank()

         ! Update the state with the new estimates
         sT_start = sT_temp
         mT_start = mT_temp
         !if (jacobian) then
         !ds_start = ds_temp
         !dm_start = dm_temp
         !end if
         call f_debug('sT_start                          ', sT_start)

         ! Aggregate data from substep to timestep
         ! Record the new values of ST
         STcum_topbot_start(:, 0) = STcum_topbot_start(:, 1)
         do concurrent(c=0:total_num_substeps - 1)
            jt_c = jt_substep_at_(c)
            if (jt_c < total_num_substeps) then
               STcum_topbot_start(jt_c + 1, 1) = STcum_topbot_start(jt_c + 1, 0) + sT_start(c)*dt_substep
            end if
         end do
         STcum_topbot_start(0, 1) = STcum_topbot_start(0, 0) + sT_init_fullstep(iT_fullstep)*dt_substep
         call f_debug('STcum_topbot_start t              ', STcum_topbot_start(:, 0))
         call f_debug('STcum_topbot_start b              ', STcum_topbot_start(:, 1))

         ! update output conc and old water frac
         do concurrent(jt_fullstep=0:timeseries_length - 1, jt_is_which_substep=0:n_substeps - 1, &
                       iq=0:numflux - 1, Q_fullstep(jt_fullstep, iq) .gt. 0)
            c = mod(total_num_substeps + jt_fullstep*n_substeps + jt_is_which_substep - iT_substep, total_num_substeps)
            C_Q_fullstep(jt_fullstep, iq, :) = &
               C_Q_fullstep(jt_fullstep, iq, :) + (mQ_aver(c, iq, :)*dt_substep/Q_fullstep(jt_fullstep, iq))/n_substeps
         end do
         do concurrent(c=0:total_num_substeps - 1)
            jt_c = jt_fullstep_at_(c)
            P_old_fullstep(jt_c, :) = P_old_fullstep(jt_c, :) - pQ_aver(c, :)*dt_substep/n_substeps
         end do

         ! Get the timestep-averaged transit time distribution
         if (iT_fullstep < max_age - 1) then
            do concurrent(outputstep=0:num_output_fullsteps - 1, jt_is_which_substep=0:n_substeps - 1, &
                          jt_is_which_substep < substep)
               jt_fullstep = output_these_fullsteps(outputstep)
               c = mod(total_num_substeps + jt_fullstep*n_substeps + jt_is_which_substep - iT_substep, total_num_substeps)
               pQ_outputstep(outputstep, :, iT_fullstep + 1) = &
                  pQ_outputstep(outputstep, :, iT_fullstep + 1) + pQ_aver(c, :)*norm
               mQ_outputstep(outputstep, :, :, iT_fullstep + 1) = &
                  mQ_outputstep(outputstep, :, :, iT_fullstep + 1) + mQ_aver(c, :, :)*norm
               mR_outputstep(outputstep, :, iT_fullstep + 1) = &
                  mR_outputstep(outputstep, :, iT_fullstep + 1) + mR_aver(c, :)*norm
            end do
         end if

         do concurrent(outputstep=0:num_output_fullsteps - 1, jt_is_which_substep=0:n_substeps - 1, &
                       jt_is_which_substep .ge. substep)
            jt_fullstep = output_these_fullsteps(outputstep)
            c = mod(total_num_substeps + jt_fullstep*n_substeps + jt_is_which_substep - iT_substep, total_num_substeps)
            pQ_outputstep(outputstep, :, iT_fullstep) = &
               pQ_outputstep(outputstep, :, iT_fullstep) + pQ_aver(c, :)*norm
            mQ_outputstep(outputstep, :, :, iT_fullstep) = &
               mQ_outputstep(outputstep, :, :, iT_fullstep) + mQ_aver(c, :, :)*norm
            mR_outputstep(outputstep, :, iT_fullstep) = &
               mR_outputstep(outputstep, :, iT_fullstep) + mR_aver(c, :)*norm
         end do

         !if (jacobian) then
         !do outputstep = 0, num_output_fullsteps - 1
         !jt_fullstep = output_these_fullsteps(outputstep)
         !do jt_substep = 0, n_substeps - 1
         !c = mod(total_num_substeps + jt_fullstep*n_substeps + jt_substep - iT_substep, total_num_substeps)
         !do iq = 0, numflux - 1
         !if (Q_fullstep(jt_fullstep, iq) > 0) then
         !dW_outputstep(outputstep, :, iq) = dW_outputstep(outputstep, :, iq) + fsQ_aver(c, :, iq)/Q_fullstep(jt_fullstep, iq)*norm*dt
         !do ic = component_index_list(iq), component_index_list(iq + 1) - 1
         !do ip = args_index_list(ic), args_index_list(ic + 1) - 1
         !dW_outputstep(outputstep, ip, iq) = dW_outputstep(outputstep, ip, iq) + fs_aver(c, ip)/Q_fullstep(jt_fullstep, iq)*norm*dt
         !end do
         !end do
         !dC_outputstep(outputstep, :, iq, :) = dC_outputstep(outputstep, :, iq, :) &
         !+ fmQ_aver(c, :, iq, :)/Q_fullstep(jt_fullstep, iq)*norm*dt
         !do ic = component_index_list(iq), component_index_list(iq + 1) - 1
         !do ip = args_index_list(ic), args_index_list(ic + 1) - 1
         !dC_outputstep(outputstep, ip, iq, :) = dC_outputstep(outputstep, ip, iq, :) &
         !+ fm_aver(c, ip, :)/Q_fullstep(jt_fullstep, iq)*norm*dt
         !end do
         !end do
         !end if
         !end do
         !end do
         !end do
         !end if

         do concurrent(outputstep=0:num_output_fullsteps - 1)
            jt_fullstep = output_these_fullsteps(outputstep)
            ! Extract substep state at timesteps
            ! age-ranked storage at the end of the timestep
            jt_is_which_substep = n_substeps - 1
            c = mod(total_num_substeps + jt_fullstep*n_substeps + jt_is_which_substep - iT_substep, total_num_substeps)
            sT_outputstep(outputstep + 1, iT_fullstep) = &
               sT_outputstep(outputstep + 1, iT_fullstep) + sT_start(c)/n_substeps
            ! parameter sensitivity
            ! Age-ranked solute mass
            mT_outputstep(outputstep + 1, :, iT_fullstep) = &
               mT_outputstep(outputstep + 1, :, iT_fullstep) + mT_start(c, :)/n_substeps
            ! parameter sensitivity
            !if (jacobian) then
            !ds_outputstep(outputstep + 1, :, iT_fullstep) = ds_outputstep(outputstep + 1, :, iT_fullstep) + ds_start(c, :)/n_substeps
            !dm_outputstep(outputstep + 1, :, :, iT_fullstep) = dm_outputstep(outputstep + 1, :, :, iT_fullstep) + dm_start(c, :, :)/n_substeps
            !end if
         end do

         call f_debug('sT_outputstep(iT_fullstep, :)     ', sT_outputstep(:, iT_fullstep))
         call f_debug('pQ_outputstep(iT_fullstep, :, 0)', pQ_outputstep(:, 0, iT_fullstep))

         iT_prev_fullstep = iT_fullstep

      end do
      ! End of the substep loop

      ! Print some updates
      if (mod(iT_fullstep, 10) .eq. 0) then
         write (tempdebugstring, *) '...Done ', (iT_fullstep), &
            'of', (max_age)
         call f_verbose(tempdebugstring)
      end if
   end do
   ! End of the main agestep loop

   ! Calculate a water balance
   ! Difference of starting and ending age-ranked storage
   do iT_fullstep = 0, max_age - 1
      do outputstep = 0, num_output_fullsteps - 1
         jt_fullstep = output_these_fullsteps(outputstep)
         if (iT_fullstep == 0) then
            WaterBalance_outputstep(outputstep, iT_fullstep) = &
               J_fullstep(jt_fullstep) - sT_outputstep(outputstep + 1, iT_fullstep)
         else
            WaterBalance_outputstep(outputstep, iT_fullstep) = &
               sT_outputstep(outputstep, iT_fullstep - 1) - sT_outputstep(outputstep + 1, iT_fullstep)
         end if
         ! subtract time-averaged water fluxes
         do iq = 0, numflux - 1
            WaterBalance_outputstep(outputstep, iT_fullstep) = &
               WaterBalance_outputstep(outputstep, iT_fullstep) - &
               (Q_fullstep(jt_fullstep, iq)*pQ_outputstep(outputstep, iq, iT_fullstep))*dt
         end do

         ! Calculate a solute balance
         ! Difference of starting and ending age-ranked mass
         if (iT_fullstep == 0) then
            do s = 0, numsol - 1
               SoluteBalance_outputstep(outputstep, s, iT_fullstep) = &
                  C_J_fullstep(jt_fullstep, s)*J_fullstep(jt_fullstep) - mT_outputstep(outputstep + 1, s, iT_fullstep)*dt
            end do
         else
            SoluteBalance_outputstep(outputstep, :, iT_fullstep) = &
               mT_outputstep(outputstep, :, iT_fullstep - 1)*dt - mT_outputstep(outputstep + 1, :, iT_fullstep)*dt
         end if
         ! Subtract timestep-averaged mass fluxes
         do iq = 0, numflux - 1
            SoluteBalance_outputstep(outputstep, :, iT_fullstep) = &
               SoluteBalance_outputstep(outputstep, :, iT_fullstep) - (mQ_outputstep(outputstep, iq, :, iT_fullstep))*dt
         end do
         ! Reacted mass
         SoluteBalance_outputstep(outputstep, :, iT_fullstep) = &
            SoluteBalance_outputstep(outputstep, :, iT_fullstep) + mR_outputstep(outputstep, :, iT_fullstep)*dt
      end do

   end do ! End of main loop

   call f_verbose('...Finalizing...')

   ! From the old water concentration
   do concurrent(s=0:numsol - 1, iq=0:numflux - 1)
      where (Q_fullstep(:, iq) > 0)
         C_Q_fullstep(:, iq, s) = C_Q_fullstep(:, iq, s) + alpha_fullstep(:, iq, s)*C_old(s)*P_old_fullstep(:, iq)
      end where
   end do

   !if (jacobian) then
   !do concurrent(s=0:numsol - 1, iq=0:numflux - 1, ip=0:numargs_total - 1)
   !where (Q_fullstep(:, iq) > 0)
   !dC_outputstep(:, ip, iq, s) = dC_outputstep(:, ip, iq, s) - C_old(s)*dW_outputstep(:, ip, iq)
   !end where
   !end do
   !end if

   call f_verbose('...Finished...')

contains

   subroutine get_flux(sT, mT, pQ, mQ, mR, stepfraction)
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1) :: sT
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mT
      real(8), intent(inout), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1) :: pQ
      real(8), intent(inout), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ
      real(8), intent(inout), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mR
      real(8), intent(in) :: stepfraction

      call f_debug('GET FLUX  rk           ', (/rk*one8, iT_substep*one8/))
      call f_debug('sT                     ', sT(:))
      call f_debug('pQ      start          ', pQ(:, 0))

      call calculate_pQ(sT, pQ, stepfraction)

      ! Solute mass flux accounting
      mQ = 0.
      do concurrent(s=0:numsol - 1, iq=0:numflux - 1, c=0:total_num_substeps - 1, sT(c) .gt. 0)
         ! Get the mass flux out
         jt_c = jt_fullstep_at_(c)
         mQ(c, iq, s) = mT(c, s)*alpha_fullstep(jt_c, iq, s)*Q_fullstep(jt_c, iq) &
                        *pQ(c, iq)/sT(c)
      end do

      ! Reaction mass accounting
      ! If there are first-order reactions, get the total mass rate
      do concurrent(s=0:numsol - 1, c=0:total_num_substeps - 1, k1_fullstep(jt_fullstep_at_(c), s) .gt. 0)
         jt_c = jt_fullstep_at_(c)
         mR(c, s) = k1_fullstep(jt_c, s)*(C_eq_fullstep(jt_c, s)*sT(c) - mT(c, s))
      end do

      !if (jacobian) then
      !fs = 0.
      !fsQ = 0.
      !do concurrent(c=0:total_num_substeps - 1, iq=0:numflux - 1, sT(c) .gt. 0)
      !fsQ(c, :, iq) = fsQ(c, :, iq) &
      !+ ds(c, :)*pQ(c, iq)*Q_fullstep(jt_c, iq)/sT(c)
      !end do
      !do iq = 0, numflux - 1
      !do ic = component_index_list(iq), component_index_list(iq + 1) - 1
      !do c = 0, total_num_substeps - 1
      !jt_c = jt_fullstep_at_(c)
                  !! sensitivity to point before the start
      !if ((leftbreakpt_top(c, ic) >= 0) .and. (leftbreakpt_top(c, ic) < numargs_list(ic) - 1)) then
      !ip = args_index_list(ic) + leftbreakpt_top(c, ic)
      !call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
      !dS = SAS_args(jt_c, ip + 1) - SAS_args(jt_c, ip)
      !dP = P_list(jt_c, ip + 1) - P_list(jt_c, ip)
      !call f_debug('dP/dS start    ', (/dP/dS/))
      !fs(c, ip) = fs(c, ip) &
      !+ dP/(dS*dS)*sT(c)*weights_fullstep(jt_c, ic)*Q_fullstep(jt_c, iq)
      !end if
                  !! sensitivity to point after the end
      !if ((leftbreakpt_bot(c, ic) + 1 > 0) .and. (leftbreakpt_bot(c, ic) + 1 <= numargs_list(ic) - 1)) then
      !ip = args_index_list(ic) + leftbreakpt_bot(c, ic) + 1
      !call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
      !dS = SAS_args(jt_c, ip) - SAS_args(jt_c, ip - 1)
      !dP = P_list(jt_c, ip) - P_list(jt_c, ip - 1)
      !call f_debug('dP/dS end      ', (/dP/dS/))
      !fs(c, ip) = fs(c, ip) &
      !- dP/(dS*dS)*sT(c)*weights_fullstep(jt_c, ic)*Q_fullstep(jt_c, iq)
      !end if
                  !! sensitivity to point within
      !if (leftbreakpt_bot(c, ic) > leftbreakpt_top(c, ic)) then
      !call f_debug('leftbreakpt_bot, _start', &
      !(/leftbreakpt_bot(ic, c)*one8, leftbreakpt_top(ic, c)*one8/))
      !do leftbreakpt = leftbreakpt_top(c, ic) + 1, leftbreakpt_bot(c, ic)
      !ip = args_index_list(ic) + leftbreakpt
      !call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
      !if (leftbreakpt > 0) then
      !dSs = SAS_args(jt_c, ip) - SAS_args(jt_c, ip - 1)
      !dPs = P_list(jt_c, ip) - P_list(jt_c, ip - 1)
      !else
      !dSs = 1.
      !dPs = 0.
      !end if
      !if (leftbreakpt < numargs_list(ic) - 1) then
      !dSe = SAS_args(jt_c, ip + 1) - SAS_args(jt_c, ip)
      !dPe = P_list(jt_c, ip + 1) - P_list(jt_c, ip)
      !else
      !dSe = 1.
      !dPe = 0.
      !end if
      !call f_debug('dP/dS middle   ', (/dPe/dSe, dPs/dSs/))
      !fs(c, ip) = fs(c, ip) &
      !- (dPe/dSe - dPs/dSs)/dt_substep*weights_fullstep(jt_c, ic)*Q_fullstep(jt_c, iq)
      !end do
      !end if
      !end do
      !end do
      !end do
      !fm = 0
      !fmQ = 0
      !do s = 0, numsol - 1
      !do iq = 0, numflux - 1
      !do ip = 0, numargs_total - 1
      !do c = 0, total_num_substeps - 1
      !if (sT(c) > 0) then
      !fmQ(c, ip, iq, s) = fmQ(c, ip, iq, s) &
      !+ dm(c, ip, s)*alpha_fullstep(jt_c, iq, s)*Q_fullstep(jt_c, iq) &
      !*pQ(c, iq)/sT(c)
      !end if
      !end do
      !end do
      !end do
      !end do
      !do s = 0, numsol - 1
      !do ip = 0, numargs_total - 1
      !do c = 0, total_num_substeps - 1
      !jt_c = jt_fullstep_at_(c)
      !fmR(c, ip, s) = fmR(c, ip, s) &
      !+ k1_fullstep(jt_c, s)*(C_eq_fullstep(jt_c, s)*ds(c, ip) - dm(c, ip, s))
      !end do
      !end do
      !end do
      !do iq = 0, numflux - 1
      !do ic = component_index_list(iq), component_index_list(iq + 1) - 1
      !do c = 0, total_num_substeps - 1
      !jt_c = jt_fullstep_at_(c)
                  !! sensitivity to point before the start
      !if ((leftbreakpt_top(c, ic) >= 0) .and. (leftbreakpt_top(c, ic) < numargs_list(ic) - 1)) then
      !ip = args_index_list(ic) + leftbreakpt_top(c, ic)
      !dS = SAS_args(jt_c, ip + 1) - SAS_args(jt_c, ip)
      !dP = P_list(jt_c, ip + 1) - P_list(jt_c, ip)
      !fm(c, ip, :) = fm(c, ip, :) &
      !+ dP/(dS*dS)*mT(c, :) &
      !*alpha_fullstep(jt_c, iq, :)*weights_fullstep(jt_c, ic)*Q_fullstep(jt_c, iq)
      !end if
                  !! sensitivity to point after the end
      !if ((leftbreakpt_bot(c, ic) + 1 > 0) .and. (leftbreakpt_bot(c, ic) + 1 <= numargs_list(ic) - 1)) then
      !ip = args_index_list(ic) + leftbreakpt_bot(c, ic) + 1
      !dS = SAS_args(jt_c, ip) - SAS_args(jt_c, ip - 1)
      !dP = P_list(jt_c, ip) - P_list(jt_c, ip - 1)
      !fm(c, ip, :) = fm(c, ip, :) &
      !- dP/(dS*dS)*mT(c, :) &
      !*alpha_fullstep(jt_c, iq, :)*weights_fullstep(jt_c, ic)*Q_fullstep(jt_c, iq)
      !end if
                  !! sensitivity to point within
      !if (leftbreakpt_bot(c, ic) > leftbreakpt_top(c, ic)) then
      !do leftbreakpt = leftbreakpt_top(c, ic) + 1, leftbreakpt_bot(c, ic)
      !ip = args_index_list(ic) + leftbreakpt
      !if (leftbreakpt > 0) then
      !dSs = SAS_args(jt_c, ip) - SAS_args(jt_c, ip - 1)
      !dPs = P_list(jt_c, ip) - P_list(jt_c, ip - 1)
      !else
      !dSs = 1.
      !dPs = 0.
      !end if
      !if (leftbreakpt < numargs_list(ic) - 1) then
      !dSe = SAS_args(jt_c, ip + 1) - SAS_args(jt_c, ip)
      !dPe = P_list(jt_c, ip + 1) - P_list(jt_c, ip)
      !else
      !dSe = 1.
      !dPe = 0.
      !end if
      !fm(c, ip, :) = fm(c, ip, :) &
      !- (dPe/dSe - dPs/dSs)*mT(c, :)/sT(c)/dt_substep &
      !*weights_fullstep(jt_c, ic)*Q_fullstep(jt_c, iq)
      !end do
      !end if
      !end do
      !end do
      !end do
      !end if
      call f_debug('pQ_temp end            ', pQ_temp(:, 0))
   end subroutine get_flux

   subroutine add_to_average(pQ, mQ, mR, coeff)
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1) :: pQ
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mR
      real(8), intent(in) :: coeff
      ! Average the rates of change using weights according to Runge-Kutta algorithm
      pQ_aver = pQ_aver + coeff*pQ
      mQ_aver = mQ_aver + coeff*mQ
      mR_aver = mR_aver + coeff*mR
      !if (jacobian) then
      !fs_aver = fs_aver   + coeff*fs
      !fsQ_aver = fsQ_aver + coeff*fsQ
      !fm_aver = fm_aver   + coeff*fm
      !fmR_aver = fmR_aver + coeff*fmR
      !fmQ_aver = fmQ_aver + coeff*fmQ
      !end if
   end subroutine add_to_average

   subroutine new_state(sT, mT, pQ, mQ, mR, stepfraction)
      real(8), intent(inout), dimension(0:timeseries_length*n_substeps - 1) :: sT
      real(8), intent(inout), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mT
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1) :: pQ
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1, 0:numsol - 1) :: mR
      real(8), intent(in) :: stepfraction

      call f_debug('NEW STATE rk           ', (/rk*one8, iT_substep*one8/))
      call f_debug('pQ                ', pQ(:, 0))
      dt_numerical_solution = dt_substep*stepfraction

      ! Calculate the new age-ranked storage
      sT = sT_start ! Initial value
      mT = mT_start + mR*dt_numerical_solution ! Initial value + reaction
      call f_debug('sT 1              ', sT(:))

      ! Fluxes in & out
      if (iT_substep == 0) then
         do concurrent(c=0:total_num_substeps - 1)
            jt_c = jt_fullstep_at_(c)
            sT(c) = sT(c) + J_fullstep(jt_c)*stepfraction
            mT(c, :) = mT(c, :) + J_fullstep(jt_c)*C_J_fullstep(jt_c, :)*stepfraction
         end do
         call f_debug('sT 2              ', sT(:))
      end if

      do concurrent(c=0:total_num_substeps - 1)
         jt_c = jt_fullstep_at_(c)
         sT(c) = sT(c) - sum(Q_fullstep(jt_c, :)*pQ(c, :))*dt_numerical_solution
         if (sT(c) < 0) then
            sT(c) = 0
         end if
      end do

      mT = mT - sum(mQ, dim=2)*dt_numerical_solution

      call f_debug('sT 3              ', sT(:))

      !if (jacobian) then
         !! Calculate new parameter sensitivity
      !ds = ds_start - sum(fsQ, dim=3)*dt_numerical_solution
      !dm = dm_start - sum(fmQ, dim=3)*dt_numerical_solution
      !ds = ds - fs*dt_numerical_solution
      !dm = dm - fm*dt_numerical_solution - fmR*dt_numerical_solution
      !end if

   end subroutine new_state

   subroutine calculate_pQ(sT, pQ, stepfraction)
      implicit none
      real(8), intent(inout), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1) :: pQ
      real(8), intent(in), dimension(0:timeseries_length*n_substeps - 1) :: sT
      real(8), intent(in) :: stepfraction
      integer :: topbot
      real(8), dimension(0:timeseries_length*n_substeps - 1) :: PQcum_component
      real(8), dimension(0:timeseries_length*n_substeps - 1, 0:1) :: STcum_topbot
      real(8), dimension(0:timeseries_length*n_substeps - 1, 0:numflux - 1, 0:1) :: PQcum_topbot
      !integer :: i_, ia_, na_
      !logical :: foundit

      if ((iT_substep == 0) .and. (stepfraction == 0)) then
         STcum_topbot = 0
         PQcum_topbot = 0
         !leftbreakpt_topbot = -1
         pQ = 0
         call f_debug('pQ is zero at start   r', pQ(:, 0))
      else
         if (iT_substep == 0) then
            STcum_topbot(:, 0) = 0
         else
            do concurrent(c=0:total_num_substeps - 1)
               STcum_topbot(c, 0) = STcum_topbot_start(jt_substep_at_(c), 0)*(1 - stepfraction) &
                                    + STcum_topbot_start(jt_substep_at_(c) + 1, 1)*(stepfraction)
            end do
         end if
         STcum_topbot(:, 1) = STcum_topbot(:, 0) + sT*dt_substep

         PQcum_topbot = 0
         do topbot = 0, 1
            ! Main lookup loop
            ! This is where we calculate the SAS functions
            !do concurrent(iq=0:numflux - 1)
            !   do concurrent(ic=component_index_list(iq):component_index_list(iq + 1) - 1)
            do iq = 0, numflux - 1
               do ic = component_index_list(iq), component_index_list(iq + 1) - 1
                  PQcum_component = 0
                  select case (component_type(ic))
                  case (-1)
                     do concurrent(c=0:total_num_substeps - 1)
                        jt_c = jt_fullstep_at_(c)
                        PQcum_component(c) = &
                           piecewiselinear_SAS_function(STcum_topbot(c, topbot), &
                                                   SAS_args(args_index_list(ic):args_index_list(ic) + numargs_list(ic) - 1, jt_c), &
                                                     P_list(args_index_list(ic):args_index_list(ic) + numargs_list(ic) - 1, jt_c), &
                                               grad_precalc(args_index_list(ic):args_index_list(ic) + numargs_list(ic) - 1, jt_c), &
                                                        numargs_list(ic))
                     end do
                  case (1)
                     !Gamma distribution
                     do concurrent(c=0:total_num_substeps - 1, sT(c) > 0)
                        jt_c = jt_fullstep_at_(c)
                        PQcum_component(c) = &
                           gamma_SAS_function(STcum_topbot(c, topbot), &
                                              SAS_args(args_index_list(ic):(args_index_list(ic) + 2), jt_c))
                     end do
                  case (2)
                     !beta distribution
                     do concurrent(c=0:total_num_substeps - 1, sT(c) > 0)
                        jt_c = jt_fullstep_at_(c)
                        PQcum_component(c) = &
                           beta_SAS_function(STcum_topbot(c, topbot), &
                                             SAS_args(args_index_list(ic):(args_index_list(ic) + 3), jt_c))
                     end do
                  case (3)
                     !kumaraswamy distribution
                     do concurrent(c=0:total_num_substeps - 1, sT(c) > 0)
                        jt_c = jt_fullstep_at_(c)
                        PQcum_component(c) = &
                           kumaraswamy_SAS_function(STcum_topbot(c, topbot), &
                                                    SAS_args(args_index_list(ic):(args_index_list(ic) + 3), jt_c))
                     end do
                  end select
                  call f_debug('STcum_topbot        ', STcum_topbot(:, topbot))
                  call f_debug('PQcum_component        ', PQcum_component)
                  do concurrent(c=0:total_num_substeps - 1, sT(c) > 0)
                     jt_c = jt_fullstep_at_(c)
                     PQcum_topbot(c, iq, topbot) = &
                        PQcum_topbot(c, iq, topbot) + weights_fullstep(jt_c, ic)*PQcum_component(c)
                  end do
               end do
            end do
         end do

         pQ = (PQcum_topbot(:, :, 1) - PQcum_topbot(:, :, 0))/dt_substep
         call f_debug('pQ is not zero, right? ', pQ(:, 0))

      end if

      do iq = 0, numflux - 1
         where (sT == 0)
            pQ(:, iq) = 0
         end where
      end do

   end subroutine calculate_pQ

   real(8) PURE FUNCTION piecewiselinear_SAS_function(ST, ST_breakpt, P_breakpt, grad, na_)
      real(8), INTENT(IN) :: ST
      integer, INTENT(IN) :: na_
      real(8), INTENT(IN), dimension(0:na_ - 1) :: ST_breakpt, P_breakpt, grad
      integer :: i_, ia_
      logical :: foundit
      if ((ST) .le. ST_breakpt(0)) then
         piecewiselinear_SAS_function = P_breakpt(0)
         !leftbreakpt_topbot(c, ic, topbot) = -1
      else
         ia_ = 0
         foundit = .FALSE.
         do i_ = 0, na_ - 1
            if (ST .lt. ST_breakpt(i_)) then
               ia_ = i_ - 1
               foundit = .TRUE.
               exit
            end if
         end do
         if (.not. foundit) then
            ia_ = na_ - 1
            piecewiselinear_SAS_function = P_breakpt(ia_)
         else
            piecewiselinear_SAS_function = P_breakpt(ia_) + (ST - ST_breakpt(ia_))*grad(ia_)
         end if
         !leftbreakpt_topbot(c, ic, topbot) = ia_
      end if
   end FUNCTION piecewiselinear_SAS_function

   real(8) PURE FUNCTION kumaraswamy_SAS_function(ST, params)
      real(8), INTENT(IN) :: ST
      real(8), INTENT(IN), dimension(4) :: params
      real(8) :: loc_, scale_, a_arg, b_arg, X
      loc_ = params(1)
      scale_ = params(2)
      a_arg = params(3)
      b_arg = params(4)
      X = (ST - loc_)/scale_
      X = MIN(MAX(0.0, X), 1.0)
      kumaraswamy_SAS_function = 1 - (1 - X**a_arg)**b_arg
   end function kumaraswamy_SAS_function

   real(8) PURE FUNCTION beta_SAS_function(ST, params)
      real(8), INTENT(IN) :: ST
      real(8), INTENT(IN), dimension(4) :: params
      real(8) :: loc_, scale_, a_arg, b_arg, X
      loc_ = params(1)
      scale_ = params(2)
      a_arg = params(3)
      b_arg = params(4)
      X = (ST - loc_)/scale_
      X = MIN(MAX(0.0, X), 1.0)
      beta_SAS_function = cum_beta_pure(X, a_arg, b_arg)
   end function beta_SAS_function

   real(8) PURE FUNCTION gamma_SAS_function(ST, params)
      real(8), INTENT(IN) :: ST
      real(8), INTENT(IN), dimension(3) :: params
      real(8) :: loc_, scale_, a_arg, X
      loc_ = params(1)
      scale_ = params(2)
      a_arg = params(3)
      X = (ST - loc_)/scale_
      gamma_SAS_function = cum_gamma_pure(X, a_arg)
   end function gamma_SAS_function

   subroutine f_debug_blank()
      ! Prints a blank line
      if (debug) then
         print *, ''
      end if
   end subroutine f_debug_blank

   subroutine f_debug(debugstring, debugdblepr)
      ! Prints debugging information
      implicit none
      character(len=*), intent(in) :: debugstring
      real(8), dimension(:), intent(in) :: debugdblepr
      if (debug) then
         print 1, debugstring, debugdblepr
1        format(A26, *(f16.10))
      end if
   end subroutine f_debug

   subroutine f_warning(debugstring)
      ! Prints informative information
      implicit none
      character(len=*), intent(in) :: debugstring
      if (warning) then
         print *, debugstring
      end if
   end subroutine f_warning

   subroutine f_verbose(debugstring)
      ! Prints informative information
      implicit none
      character(len=*), intent(in) :: debugstring
      if (verbose) then
         print *, debugstring
      end if
   end subroutine f_verbose

end subroutine solveSAS
