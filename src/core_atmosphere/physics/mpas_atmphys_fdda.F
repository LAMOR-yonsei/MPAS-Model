! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_atmphys_fdda

    use mpas_derived_types, only : mpas_pool_type, mpas_clock_type, block_type, mpas_time_type, mpas_timeInterval_type, MPAS_NOW, &
                                   MPAS_STREAM_LATEST_BEFORE, MPAS_STREAM_EARLIEST_STRICTLY_AFTER, MPAS_STREAM_NEAREST, &
                                   MPAS_streamManager_type
    use mpas_pool_routines, only : mpas_pool_get_array, mpas_pool_get_dimension, mpas_pool_get_subpool, mpas_pool_shift_time_levels, mpas_pool_get_config
    use mpas_kind_types, only : RKIND, StrKIND
    use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_timeInterval, mpas_set_time, operator(-)
    use mpas_timer, only : mpas_timer_start, mpas_timer_stop

    use mpas_atmphys_constants
    use mpas_atmphys_utilities
    use mpas_atmphys_vars

    public :: mpas_atm_get_fdda_tend, &
              mpas_atm_apply_fdda
    private
    type (MPAS_Time_Type) :: FDDA_intv_end
    integer, private:: i,k

    contains

    !***********************************************************************
    !
    !  routine mpas_atm_get_fdda_tend
    !
    !> \brief   Reads variables for FDDA
    !> \author  Sang-Hun Park
    !> \date    22 August 2020
    !> \details 
    !>  This routine reads from the 'fdda_in' stream all variables in the 'fdda'
    !>  pool. In this time, coupled momentun (ru/rv), potential temperature
    !>  (rtheta) and qvapor (rqv) are available. 
    !
    !-----------------------------------------------------------------------

   subroutine mpas_atm_get_fdda_tend(clock, streamManager, block, firstCall, ierr)

        use mpas_stream_manager, only : mpas_stream_mgr_read
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only : MPAS_STREAM_MGR_NOERR, MPAS_LOG_ERR
        use mpas_timekeeping, only : mpas_get_time
        use mpas_constants, only : rvord

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (MPAS_streamManager_type), intent(inout) :: streamManager
        type (block_type), intent(inout) :: block
        logical, intent(in) :: firstCall
        integer, intent(out) :: ierr

        character(len=StrKIND) :: fdda_intv_start_string
        character(len=StrKIND) :: fdda_intv_end_string

        type (mpas_pool_type), pointer :: mesh
        type (mpas_pool_type), pointer :: fdda_state

        real (kind=RKIND) :: dt

        real (kind=RKIND), dimension(:,:), pointer :: u, v, theta
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars

        real (kind=RKIND), dimension(:,:), pointer :: fdda_tend_u, fdda_tend_v, fdda_tend_theta
        real (kind=RKIND), dimension(:,:,:), pointer :: fdda_tend_scalars
        integer, dimension(:,:), pointer :: cellsOnEdge
        integer, pointer  :: nCells, nEdges, nVertLevels, index_qv

        integer :: dd_intv, s_intv, sn_intv, sd_intv
        type (MPAS_Time_Type) :: currTime
        type (MPAS_TimeInterval_Type) :: fdda_interval
        character(len=StrKIND) :: read_time

        integer :: iEdge
        integer :: cell1, cell2

        ierr = 0

        call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
        call mpas_pool_get_subpool(block % structs, 'fdda', fdda_state)

        if (firstCall) then
            call MPAS_stream_mgr_read(streamManager, streamID='fdda_in', timeLevel=2, whence=MPAS_STREAM_LATEST_BEFORE, &
                                      actualWhen=read_time, ierr=ierr)
            if (ierr /= MPAS_STREAM_MGR_NOERR) then
               call mpas_log_write('Could not read from ''fdda_in'' stream on or before the current date '// &
                                   'to update lateral boundary tendencies', messageType=MPAS_LOG_ERR)
               ierr = 1
            end if
        else
            call mpas_pool_shift_time_levels(fdda_state)
            call MPAS_stream_mgr_read(streamManager, streamID='fdda_in', timeLevel=2, whence=MPAS_STREAM_EARLIEST_STRICTLY_AFTER, &
                                      actualWhen=read_time, ierr=ierr)
            if (ierr /= MPAS_STREAM_MGR_NOERR) then
               call mpas_log_write('Could not read from ''fdda_in'' stream after the current date '// &
                                   'to update lateral boundary tendencies', messageType=MPAS_LOG_ERR)
               ierr = 1
            end if
        end if
        if (ierr /= 0) then
           return
        end if

        call mpas_set_time(currTime, dateTimeString=trim(read_time))

        !call MPAS_stream_mgr_read(streamManager, streamID='fdda_in',
        !whence=MPAS_STREAM_NEAREST, actualWhen=read_time, ierr=ierr)
        if (ierr /= MPAS_STREAM_MGR_NOERR) then
           call mpas_log_write('Could not read from ''fdda_in'' stream on or before the current date '// &
                               'to save variables for FDDA', messageType=MPAS_LOG_ERR)
           ierr = 1
        end if

        if (ierr /= 0) then
           return
        end if

        call mpas_pool_get_array(fdda_state, 'fdda_u', u, 2)
        call mpas_pool_get_array(fdda_state, 'fdda_v', v, 2)
        call mpas_pool_get_array(fdda_state, 'fdda_theta', theta, 2)
        call mpas_pool_get_array(fdda_state, 'fdda_scalars', scalars, 2)

        call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(fdda_state, 'index_qv', index_qv)

        if (.not. firstCall) then
            fdda_interval = currTime - FDDA_intv_end
            call mpas_get_timeInterval(interval=fdda_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
            dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) &
                 + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))

            call mpas_pool_get_array(fdda_state, 'fdda_u', u, 2)
            call mpas_pool_get_array(fdda_state, 'fdda_v', v, 2)
            call mpas_pool_get_array(fdda_state, 'fdda_theta', theta, 2)
            call mpas_pool_get_array(fdda_state, 'fdda_scalars', scalars, 2)

            call mpas_pool_get_array(fdda_state, 'fdda_u', fdda_tend_u, 1)
            call mpas_pool_get_array(fdda_state, 'fdda_v', fdda_tend_v, 1)
            call mpas_pool_get_array(fdda_state, 'fdda_theta', fdda_tend_theta, 1)
            call mpas_pool_get_array(fdda_state, 'fdda_scalars', fdda_tend_scalars, 1)

            !dt = 1.0_RKIND / dt
            !fdda_tend_u(:,:) = (u(:,:) - fdda_tend_u(:,:)) * dt
            !fdda_tend_v(:,:) = (v(:,:) - fdda_tend_v(:,:)) * dt
            !fdda_tend_theta(:,:) = (theta(:,:) - fdda_tend_theta(:,:)) * dt
            !fdda_tend_scalars(:,:,:) = (scalars(:,:,:) - fdda_tend_scalars(:,:,:)) * dt

            ! Logging the FDDA start and end times appears to be backwards, but
            ! until the end of this function, FDDA_intv_end == the last interval
            ! time and currTime == the next interval time.
            !
            call mpas_get_time(FDDA_intv_end, dateTimeString=fdda_intv_start_string)
            call mpas_get_time(currTime, dateTimeString=fdda_intv_end_string)
            call mpas_log_write('----------------------------------------------------------------------')
            call mpas_log_write('Read variables from a FDDA file. FDDAs are now valid')
            !call mpas_log_write('Read variables from a FDDA file on
            !...'//trim(read_time))
            call mpas_log_write('from '//trim(fdda_intv_start_string)//' to '//trim(fdda_intv_end_string))
            call mpas_log_write('----------------------------------------------------------------------')

        end if

        FDDA_intv_end = currTime

    end subroutine mpas_atm_get_fdda_tend

!
! --> This routine is from O. Russell Bullock from EPA and modified including dss_fdda
!
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>      Original sourcecode by O. Russell Bullock Jr. (bullock.russell@epa.gov) / 2016-03-30. 
!>
!>    * Added various controls for length-scaled FDDA application 
!>      O. Russell Bullock Jr. (bullock.russell@epa.gov) / 2016-08-11. 
!>    * Modified to only require reads of "fdda_new" data arrays.  "fdda_old" is initially 
!>      set to initial conditions.  At end of FDDA data interval, "fdda_old" is set to 
!>      "fdda_new" and updated "fdda_new" data comes from FDDA input stream.
!>      O. Russell Bullock Jr. (bullock.russell@epa.gov) / 2016-08-30.

!==================================================================================================
 subroutine mpas_atm_apply_fdda(xtime_s,configs,mesh,state,time_lev,fdda_state,diag,diag_physics,tend_physics,its,ite)
!==================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: fdda_state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics

 integer,intent(in):: time_lev, its, ite
 real(kind=RKIND),intent(in):: xtime_s

!inout variables:
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,pointer:: index_qv
 real(kind=RKIND),pointer:: fdda_int, fdda_xnutr
 real(kind=RKIND),dimension(:,:),pointer  :: rthfddaten,rqvfddaten
 real(kind=RKIND),dimension(:,:),pointer  :: rufddaten,rvfddaten 
 real(kind=RKIND),dimension(:,:),pointer  :: qv_fdda_old,th_fdda_old,u_fdda_old,v_fdda_old
 real(kind=RKIND),dimension(:,:),pointer  :: qv_fdda_new,th_fdda_new,u_fdda_new,v_fdda_new
 real(kind=RKIND),dimension(:,:,:),pointer  :: scalars_fdda_old
 real(kind=RKIND),dimension(:,:,:),pointer  :: scalars_fdda_new
 real(kind=RKIND),dimension(:,:),pointer  :: u,v
 real(kind=RKIND),dimension(:,:),pointer  :: qv
 real(kind=RKIND),dimension(:,:),pointer  :: theta_m
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),pointer  :: dss_fdda

!local variables:
 real(kind=RKIND)  :: tfrac
 real(kind=RKIND)  :: t_target,q_target,u_target,v_target
 real(kind=RKIND)  :: sfac
 real(kind=RKIND)  :: theta_p
 real(kind=RKIND)  :: srange

 integer :: i,k

!local parameters:
 real(kind=RKIND),parameter:: rvrd   = R_v/R_d

!==================================================================================================
 write(0,*) '--- enter apply_fdda:  at xtime_s = ',xtime_s

 call mpas_pool_get_config(configs,'config_fdda_interval',fdda_int)
 call mpas_pool_get_config(configs,'config_fdda_xnutr',fdda_xnutr)

 call mpas_pool_get_dimension(state,'index_qv',index_qv)

 call mpas_pool_get_array(fdda_state,'fdda_scalars',scalars_fdda_old,1)
 qv_fdda_old => scalars_fdda_old(index_qv,:,:)
 call mpas_pool_get_array(fdda_state,'fdda_theta',th_fdda_old,1)
 call mpas_pool_get_array(fdda_state,'fdda_u', u_fdda_old,1 )
 call mpas_pool_get_array(fdda_state,'fdda_v', v_fdda_old,1 )

 call mpas_pool_get_array(fdda_state,'fdda_scalars',scalars_fdda_new,2)
 qv_fdda_new => scalars_fdda_new(index_qv,:,:)
 call mpas_pool_get_array(fdda_state,'fdda_theta',th_fdda_new,2)
 call mpas_pool_get_array(fdda_state,'fdda_u', u_fdda_new,2 )
 call mpas_pool_get_array(fdda_state,'fdda_v', v_fdda_new,2 )

 call mpas_pool_get_array(tend_physics,'rthfddaten',rthfddaten)
 call mpas_pool_get_array(tend_physics,'rqvfddaten',rqvfddaten)
 call mpas_pool_get_array(tend_physics,'rufddaten' ,rufddaten )
 call mpas_pool_get_array(tend_physics,'rvfddaten' ,rvfddaten )

 call mpas_pool_get_array(diag,'uReconstructZonal'     ,u)
 call mpas_pool_get_array(diag,'uReconstructMeridional',v)

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)

 call mpas_pool_get_array(mesh, 'dss_fdda', dss_fdda)

!calculate time fraction within FDDA target data interval
 tfrac = mod(xtime_s,fdda_int)/fdda_int

 DO i=its,ite
 DO k=kts,kte

    theta_p = theta_m(k,i) / (1._RKIND + rvrd * qv(k,i))

    t_target = (1.0-tfrac)*th_fdda_old(k,i)+tfrac*th_fdda_new(k,i)
    rthfddaten(k,i) = fdda_xnutr * dss_fdda(k,i) * ( t_target - theta_p )

    q_target = (1.0-tfrac)*qv_fdda_old(k,i)+tfrac*qv_fdda_new(k,i)
    rqvfddaten(k,i) = fdda_xnutr * dss_fdda(k,i) * ( q_target - qv(k,i) )

    u_target = (1.0-tfrac)*u_fdda_old(k,i)+tfrac*u_fdda_new(k,i)
    rufddaten(k,i) = fdda_xnutr * dss_fdda(k,i) * ( u_target - u(k,i) )

    v_target = (1.0-tfrac)*v_fdda_old(k,i)+tfrac*v_fdda_new(k,i)
    rvfddaten(k,i) = fdda_xnutr * dss_fdda(k,i) * ( v_target - v(k,i) )

 END DO
 END DO

 write(0,*) '--- exit apply_fdda'

 end subroutine mpas_atm_apply_fdda

!==================================================================================================
 end module mpas_atmphys_fdda
!==================================================================================================
