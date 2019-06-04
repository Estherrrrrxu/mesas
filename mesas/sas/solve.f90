! -*- f90 -*-
      subroutine f_solve(J_ts, Q_ts, rSAS_lookup, P_list, STcum_init_ts, dt,  &
                 verbose, debug, full_outputs,&
                 CS_init_ts, C_J_ts, alpha_ts, k1_ts, C_eq_ts, C_old, STcum_ts, PQcum_ts, &
                 WaterBalance_ts, MScum_ts, MQcum_ts, MRcum_ts, C_Q_ts, SoluteBalance_ts, &
                 n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list)
      implicit none
      logical, intent(in) :: verbose, debug, full_outputs
      integer, intent(in) ::  n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list
      real(8), intent(in) :: dt
      real(8), intent(in), dimension(0:timeseries_length-1) :: J_ts
      real(8), intent(in), dimension(0:timeseries_length-1, 0:numflux-1) :: Q_ts
      real(8), intent(in), dimension(0:nP_list-1) :: P_list
      real(8), intent(in), dimension(0:nP_list-1, 0:timeseries_length-1, 0:numflux-1) :: rSAS_lookup
      real(8), intent(in), dimension(0:max_age) :: STcum_init_ts
      real(8), intent(in), dimension(0:max_age-1,0:numsol-1) :: CS_init_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numsol-1) :: C_J_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numflux-1,0:numsol-1) &
                 :: alpha_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numsol-1) :: k1_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numsol-1) :: C_eq_ts
      real(8), intent(in), dimension(0:numsol-1) :: C_old
      real(8), intent(out), dimension(0:timeseries_length-1, 0:numflux-1, 0:numsol-1) :: C_Q_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length) :: STcum_ts
      real(8), intent(out), dimension(0:max_age-1, 0:timeseries_length-1) :: &
                 WaterBalance_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length, 0:numflux-1) :: PQcum_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length, 0:numsol-1) :: MScum_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length, 0:numflux-1, 0:numsol-1) :: MQcum_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length, 0:numsol-1) :: MRcum_ts
      real(8), intent(out), dimension(0:max_age-1, 0:timeseries_length-1, 0:numsol-1) :: SoluteBalance_ts
      integer :: k, i
      real(8) :: h
      real(8), dimension(0:max_age*n_substeps) :: STcum_init
      real(8), dimension(0:max_age*n_substeps, 0:numflux-1) :: PQcum_init
      real(8), dimension(0:max_age*n_substeps, 0:numflux-1) :: PQcum_temp
      real(8), dimension(0:max_age*n_substeps) :: STcum_temp
      real(8), dimension(0:max_age*n_substeps-1) :: sT_prev
      real(8), dimension(0:max_age*n_substeps-1) :: sT_temp
      real(8), dimension(0:max_age*n_substeps-1) :: sT_substep
      real(8), dimension(0:max_age*n_substeps) :: STcum_substep
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ1
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ2
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ3
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ4
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ_substep
      real(8), dimension(0:max_age-1, 0:numflux-1) :: pQ_ts
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ1
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ2
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ3
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ4
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ_substep
      real(8), dimension(0:max_age-1, 0:numflux-1, 0:numsol-1) :: mQ_ts
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR1
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR2
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR3
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR4
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR_substep
      real(8), dimension(0:max_age-1, 0:numsol-1) :: mR_ts
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mS_prev
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mS_substep
      real(8), dimension(0:max_age*n_substeps, 0:numsol-1) :: MScum_substep
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mS_temp
      real(8) :: one8
      character(len=128) :: debugstring
      integer is, ie, iq, s, M
      one8 = 1.0
      call f_verbose(verbose,'...Initializing arrays...')
      C_Q_ts(:, :, :) = 0.
      STcum_ts(:, :) = 0.
      WaterBalance_ts(:, :) = 0.
      PQcum_ts(:, :, :) = 0.
      MScum_ts(:, :, :) = 0.
      MQcum_ts(:, :, :, :) = 0.
      MRcum_ts(:, :, :) = 0.
      SoluteBalance_ts(:, :, :) = 0.
      STcum_init(:) = 0.
      PQcum_init(:, :) = 0.
      PQcum_temp(:,:) = 0.
      STcum_temp(:) = 0.
      sT_prev(:) = 0.
      sT_temp(:) = 0.
      sT_substep(:) = 0.
      STcum_substep(:) = 0.
      pQ1(:, :) = 0.
      pQ2(:, :) = 0.
      pQ3(:, :) = 0.
      pQ4(:, :) = 0.
      pQ_substep(:, :) = 0.
      pQ_ts(:, :) = 0.
      mQ1(:, :, :) = 0.
      mQ2(:, :, :) = 0.
      mQ3(:, :, :) = 0.
      mQ4(:, :, :) = 0.
      mQ_substep(:, :, :) = 0.
      mQ_ts(:, :, :) = 0.
      mR1(:, :) = 0.
      mR2(:, :) = 0.
      mR3(:, :) = 0.
      mR4(:, :) = 0.
      mR_substep(:, :) = 0.
      mR_ts(:, :) = 0.
      mS_prev(:, :) = 0.
      mS_substep(:, :) = 0.
      MScum_substep(:, :) = 0.
      mS_temp(:, :) = 0.
      M = max_age * n_substeps
      h = dt / n_substeps

      call f_verbose(verbose,'...Setting initial conditions...')
      do i=0,max_age-1
        is = i*n_substeps
        ie = (i+1)*n_substeps-1
        sT_substep(is:ie) = (STcum_init_ts(i+1) - STcum_init_ts(i))/n_substeps
        do s=0,numsol-1
          mS_substep(is:ie,s) = (CS_init_ts(i,s) * (STcum_init_ts(i+1) - STcum_init_ts(i)))&
                          /n_substeps
        enddo
      enddo
      STcum_init = cumsum(sT_substep, M)
      do iq=0,numflux-1
        call lookup(rSAS_lookup(:,0,iq), P_list, &
                    STcum_init, PQcum_init(:,iq), nP_list, M+1)
      enddo
      STcum_ts(1:max_age,0) = STcum_init(n_substeps:M:n_substeps)
      do iq=0,numflux-1
        PQcum_ts(:,0,iq) = PQcum_init(0:M:n_substeps, iq)
      enddo
      do s=0,numsol-1
        MScum_substep(:,s) = cumsum(mS_substep(:,s), M)
        MScum_ts(:,0,s) = MScum_substep(0:M:n_substeps, s)
      enddo
      call f_verbose(verbose,'...Starting main loop...')
      do i=0,timeseries_length-1
        pQ_ts(:, :) = 0.
        mQ_ts(:, :, :) = 0.
        mR_ts(:, :) = 0.
        do k=0,n_substeps-1
          call f_debug(debug, 'Timestep, Substep', (/ i*one8, k*one8/))
          sT_prev(1:M-1) = sT_substep(0:M-2)
          sT_prev(0) = 0.
          mS_prev(1:M-1,:) = mS_substep(0:M-2,:)
          mS_prev(0,:) = 0.
          call f_debug(debug, 'RK', (/1._8/))
          call get_flux(i, sT_prev, STcum_temp, pQ1, PQcum_temp, J_ts, Q_ts, rSAS_lookup, &
                 P_list, mS_prev, mQ1, mR1, alpha_ts, k1_ts, C_eq_ts, n_substeps, &
                 numflux, numsol, max_age, timeseries_length, nP_list)
          call new_state(i, h/2, sT_prev, sT_temp, pQ1, J_ts, Q_ts, mS_prev, mS_temp, mQ1, &
                 mR1, C_J_ts, n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list)
          call f_debug(debug, 'RK', (/2._8/))
          call get_flux(i, sT_temp, STcum_temp, pQ2, PQcum_temp, J_ts, Q_ts, rSAS_lookup, &
                 P_list, mS_temp, mQ2, mR2, alpha_ts, k1_ts, C_eq_ts, n_substeps, &
                 numflux, numsol, max_age, timeseries_length, nP_list)
          call new_state(i, h/2, sT_prev, sT_temp, pQ2, J_ts, Q_ts, mS_prev, mS_temp, mQ2, &
                 mR2, C_J_ts, n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list)
          call f_debug(debug, 'RK', (/3._8/))
          call get_flux(i, sT_temp, STcum_temp, pQ3, PQcum_temp, J_ts, Q_ts, rSAS_lookup, &
                 P_list, mS_temp, mQ3, mR3, alpha_ts, k1_ts, C_eq_ts, n_substeps, &
                 numflux, numsol, max_age, timeseries_length, nP_list)
          call new_state(i, h, sT_prev, sT_temp, pQ3, J_ts, Q_ts, mS_prev, mS_temp, mQ3, mR3&
                 , C_J_ts, n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list)
          call f_debug(debug, 'RK', (/4._8/))
          call get_flux(i, sT_temp, STcum_temp, pQ4, PQcum_temp, J_ts, Q_ts, rSAS_lookup, &
                 P_list, mS_temp, mQ4, mR4, alpha_ts, k1_ts, C_eq_ts, n_substeps, &
                 numflux, numsol, max_age, timeseries_length, nP_list)
          pQ_substep = (pQ1 + 2*pQ2 + 2*pQ3 + pQ4) / 6.
          mQ_substep = (mQ1 + 2*mQ2 + 2*mQ3 + mQ4) / 6.
          mR_substep = (mR1 + 2*mR2 + 2*mR3 + mR4) / 6.
          do iq=0,numflux-1
            if (Q_ts(i,iq)==0) then
              pQ_substep(:,iq) = 0.
            endif
          enddo
          call new_state(i, h, sT_prev, sT_substep, pQ_substep, J_ts, Q_ts, mS_prev, mS_substep, mQ_substep, mR_substep&
                 , C_J_ts, n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list)
          do iq=0,numflux-1
            call f_debug(debug, 'pQ_substep', pQ_substep(:,iq))
          enddo
          do s=0,numsol-1
            call f_debug(debug, 'mS_substep', mS_substep(:,s))
            call f_debug(debug, 'mR_substep', mR_substep(:,s))
            do iq=0,numflux-1
              call f_debug(debug, 'mQ_substep', mQ_substep(:,iq,s))
            enddo
          enddo
          do iq=0,numflux-1
            do s=0,numsol-1
              if (Q_ts(i,iq)>0) then
                C_Q_ts(i,iq,s) = C_Q_ts(i,iq,s) + sum(mQ_substep(:,iq,s)) / &
                 (dt * Q_ts(i,iq)) / n_substeps
              call f_debug(debug, 'C_Q_ts', C_Q_ts(i:i,iq,s))
              endif
            enddo
          enddo
          if (full_outputs) then
            do iq=0,numflux-1
              pQ_ts(0, iq) = pQ_ts(0, iq) &
                           + sum(pQ_substep(0:k, iq))/n_substeps
              pQ_ts(1:max_age-1, iq) = pQ_ts(1:max_age-1, iq) &
                   + sum(reshape(pQ_substep(k+1:M-n_substeps+k, iq),&
                                 [n_substeps,max_age-1]),1)/n_substeps
              do s=0,numsol-1
                mQ_ts(0, iq, s) = mQ_ts(0, iq, s) &
                           + sum(mQ_substep(0:k, iq, s))/n_substeps
                mQ_ts(1:max_age-1, iq, s) = mQ_ts(1:max_age-1, iq, s)&
                           + sum(reshape(mQ_substep(k+1:M-n_substeps+k, iq, s), &
                                 [n_substeps,max_age-1]),1)/n_substeps
              enddo
            enddo
            do s=0,numsol-1
              mR_ts(0, s) = mR_ts(0, s) &
                           + sum(mR_substep(0:k, s))/n_substeps
              mR_ts(1:max_age-1, s) = mR_ts(1:max_age-1, s) &
                           + sum(reshape(mR_substep(k+1:M-n_substeps+k, s), &
                                 [n_substeps,max_age-1]),1)/n_substeps
            enddo
          endif
        enddo
        if (full_outputs) then
          STcum_substep = cumsum(sT_substep, M)
          STcum_ts(1:max_age, i+1) = STcum_substep(n_substeps:M:n_substeps)
          WaterBalance_ts(1:max_age-1, i) = diff(STcum_ts(0:max_age-1, i), max_age) - &
                   diff(STcum_ts(1:max_age, i+1), max_age)
          WaterBalance_ts(0, i) = J_ts(i) * dt - STcum_ts(1, i+1)
          do iq=0,numflux-1
            PQcum_ts(0:max_age,i+1,iq) = cumsum(pQ_ts(0:max_age-1,iq), max_age)
            WaterBalance_ts(1:max_age-1, i) = WaterBalance_ts(1:max_age-1, i) - dt &
                   * (Q_ts(i,iq) * diff(PQcum_ts(1:max_age,i+1,iq), max_age))
            WaterBalance_ts(0, i) = WaterBalance_ts(0, i) - dt * Q_ts(i,iq) * &
                   (PQcum_ts(1,i+1,iq) - PQcum_ts(0,i+1,iq))
          enddo
          do s=0,numsol-1
            MScum_substep(:,s) = cumsum(mS_substep(:,s), M)
            MScum_ts(:,i+1,s) = MScum_substep(0:M:n_substeps,s)
            MRcum_ts(:,i+1,s) = cumsum(mR_ts(:,s), max_age)
            do iq=0,numflux-1
              MQcum_ts(:,i+1,iq,s) = cumsum(mQ_ts(:,iq,s),&
                   max_age)
              C_Q_ts(i,iq,s) = C_Q_ts(i,iq,s) + alpha_ts(i,iq,s) * C_old(s) * (1 &
                   - PQcum_ts(max_age, i+1, iq))
            enddo
          enddo
          do s=0,numsol-1
            SoluteBalance_ts(1:max_age-1,i,s)=(&
                   diff(MScum_ts(0:max_age-1,i,s), max_age) - &
                   diff(MScum_ts(1:max_age,i+1,s), max_age) &
                   + dt * diff(MRcum_ts(1:max_age-1,i+1,s), max_age))
            SoluteBalance_ts(0,i,s) = (C_J_ts(i,s) * J_ts(i) * dt - MScum_ts(0,i+1,s) + &
                   dt * (MRcum_ts(0,i+1,s) - MRcum_ts(0,i+1,s)))
            do iq=0,numflux-1
              SoluteBalance_ts(0:max_age-1,i,s) = (SoluteBalance_ts(0:max_age-1,i,s&
                   ) - dt * diff(MQcum_ts(0:max_age-1,i+1,iq,s), max_age+1))
              SoluteBalance_ts(0,i,s) = (SoluteBalance_ts(0,i,s) - dt * (MQcum_ts(0,&
                   i +1,iq,s) - MQcum_ts(0,i+1,iq,s)))
            enddo
          enddo
        endif
        if (mod(i,1000).eq.1000) then
            write (debugstring,*) '...Done ', char(i), &
                                  'of', char(timeseries_length)
            call f_verbose(verbose, debugstring)
        endif
      enddo
      call f_verbose(verbose,'...Finished...')
      contains

      subroutine f_debug(debug, debugstring, debugdblepr)
      implicit none
      logical, intent(in) :: debug
      character(len=*), intent(in) :: debugstring
      real(8), dimension(:), intent(in) :: debugdblepr
      if (debug) then
          print *, debugstring, debugdblepr
          !if (size(debugdblepr)==0) then
          !    call dblepr(debugstring, -1, 0, 0)
          !else
          !    call dblepr(debugstring, -1, debugdblepr, size(debugdblepr))
          !endif
      endif
      end subroutine f_debug

      subroutine f_verbose(verbose, debugstring)
      implicit none
      logical, intent(in) :: verbose
      character(len=*), intent(in) :: debugstring
      if (verbose) then
          print *, debugstring
          !call dblepr(debugstring, -1, 0, 0)
      endif
      end subroutine f_verbose

      subroutine lookup_forwarddifference(xa, ya, x, y, dx, dy, na, n, i)
      implicit none
      integer, intent(in) ::  na, n, i
      real(8), intent(in), dimension(0:na-1) :: xa
      real(8), intent(in), dimension(0:na-1) :: ya
      real(8), intent(in), dimension(0:n-1) :: x
      real(8), intent(in), dimension(0:n-1) :: y
      real(8), intent(in), dimension(0:n-1) :: dx
      real(8), intent(out), dimension(0:n-1) :: dy
      integer :: k, j, k0
      real(8) :: dif, grad
      logical :: foundit
      k0 = 0
      do j=0,n-1
        if (x(j).le.xa(0)) then
            y(j) = ya(0)
        else
            foundit = .FALSE.
            do k=k0,na-1
                if (x(j).lt.xa(k)) then
                    k0 = k-1
                    foundit = .TRUE.
                    exit
                endif
            enddo
            if (.not. foundit) then
                y(j) = ya(na-1)
            else
                dif = x(j) - xa(k0)
                grad = (ya(k0+1)-ya(k0))/(xa(k0+1)-xa(k0))
                y(j) = ya(k0) + dif * grad
            endif
        endif
      enddo
      end subroutine

      function diff(arr, n)
      implicit none
      integer, intent(in) ::  n
      real(8), intent(in), dimension(0:n-1) :: arr
      real(8), dimension(0:n-1-1) :: diff
      integer :: i
      do i=0,n-1-1
        diff(i) = arr(i+1) - arr(i)
      enddo
      end function diff

      function cumsum(arr, n)
      implicit none
      integer, intent(in) ::  n
      real(8), intent(in), dimension(0:n-1) :: arr
      real(8), dimension(0:n) :: cumsum
      integer :: i
      cumsum(1) = 0
      do i=0,n-1
        cumsum(i+1) = cumsum(i) + arr(i)
      enddo
      end function cumsum

      subroutine get_flux(i, sT_temp, STcum_temp, pQ_temp, PQcum_temp, J_ts, Q_ts, rSAS_lookup&
                 , P_list, mS_temp, mQ_temp, mR_temp, alpha_ts, k1_ts, C_eq_ts, n_substeps, &
                 numflux, numsol, max_age, timeseries_length, nP_list)
      implicit none
      integer, intent(in) :: i
      integer, intent(in) :: n_substeps, numflux, numsol, max_age, timeseries_length, nP_list
      real(8), intent(in), dimension(0:timeseries_length-1) :: J_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numflux-1) :: Q_ts
      real(8), intent(in), dimension(0:nP_list-1) :: P_list
      real(8), intent(in), dimension(0:nP_list-1, 0:timeseries_length-1, 0:numflux-1) :: rSAS_lookup
      real(8), intent(in), dimension(0:timeseries_length-1,0:numflux-1,0:numsol-1) :: alpha_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numsol-1) :: k1_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numsol-1) :: C_eq_ts
      real(8), intent(inout), dimension(0:timeseries_length*n_substeps, 0:numflux-1) :: PQcum_temp
      real(8), intent(inout), dimension(0:timeseries_length*n_substeps) :: STcum_temp
      real(8), intent(out), dimension(0:timeseries_length*n_substeps-1, 0:numflux-1) :: pQ_temp
      real(8), intent(in), dimension(0:timeseries_length*n_substeps-1) :: sT_temp
      real(8), intent(out), dimension(0:timeseries_length*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ_temp
      real(8), intent(out), dimension(0:timeseries_length*n_substeps-1, 0:numsol-1) :: mR_temp
      real(8), intent(in), dimension(0:timeseries_length*n_substeps-1, 0:numsol-1) :: mS_temp
      integer M, iq, s
      M = max_age * n_substeps
      N = timeseries_length * n_substeps
      call f_debug(debug, 'sT_temp', sT_temp)
      call f_debug(debug, 'STcum_temp', STcum_temp)
      do s=0,numsol-1
        call f_debug(debug, 'mS_temp', mS_temp(:,s))
      enddo
      do iq=0,numflux-1
        call lookup(rSAS_lookup(:,i,iq), P_list, STcum_temp, &
                 sT_temp, pQ_temp(:,iq), nP_list, N, i)
        call f_debug(debug, 'PQcum_temp', PQcum_temp(:,iq))
        call f_debug(debug, 'pQ_temp', pQ_temp(:,iq))
        do s=0,numsol-1
          where (sT_temp(:)>0)
            mQ_temp(:,iq,s) = mS_temp(:,s) / sT_temp(:) * alpha_ts(i,iq,s) &
                 * Q_ts(i,iq) * pQ_temp(:,iq)
          elsewhere
            mQ_temp(:,iq,s) = 0.
          end where
        enddo
      enddo
      do s=0,numsol-1
         mR_temp(:,s) = k1_ts(i,s) * (C_eq_ts(i,s) * sT_temp - mS_temp(:,s))
         call f_debug(debug, 'mEt', mQ_temp(:,1,s))
         call f_debug(debug, 'mQ_temp', mQ_temp(:,0,s))
         call f_debug(debug, 'mR_temp', mR_temp(:,s))
      enddo
      end subroutine get_flux

      subroutine new_state(i, hr, sT_prev, sT_temp, pQ_temp, J_ts, Q_ts, mS_prev, mS_temp, mQ_temp, &
                 mR_temp, C_J_ts, n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list)
      integer, intent(in) :: n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list
      integer, intent(in) :: i
      real(8), intent(in) :: hr
      real(8), intent(in), dimension(0:timeseries_length-1) :: J_ts
      real(8), intent(in), dimension(0:timeseries_length-1, 0:numflux-1) :: Q_ts
      real(8), intent(inout), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ_temp
      real(8), intent(inout), dimension(0:max_age*n_substeps-1) :: sT_temp
      real(8), intent(inout), dimension(0:max_age*n_substeps-1) :: sT_prev
      real(8), intent(inout), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ_temp
      real(8), intent(inout), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mS_temp
      real(8), intent(inout), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mS_prev
      real(8), intent(inout), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR_temp
      real(8), intent(in), dimension(0:timeseries_length-1, 0:numsol-1) :: C_J_ts
      ! This needs to be done here: STcum_temp = cumsum(sT_temp, M)
      M = max_age * n_substeps
      sT_temp = sT_prev
      do iq=0,numflux-1
        sT_temp = sT_temp - Q_ts(i,iq) * pQ_temp(:,iq) * hr
      enddo
      sT_temp(0) = sT_temp(0) + J_ts(i) * hr
      do s=0,numsol-1
        mS_temp(:,s) = mS_prev(:,s) + mR_temp(:,s) * hr
        do iq=0,numflux-1
          mS_temp(:,s) = mS_temp(:,s) - mQ_temp(:,iq,s) * hr
        enddo
        mS_temp(0,s) = mS_temp(0,s) + J_ts(i) * C_J_ts(i,s) * hr
      enddo
      end subroutine new_state

      end subroutine f_solve
