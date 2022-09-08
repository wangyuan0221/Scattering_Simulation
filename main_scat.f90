! -------fortran program------- !
! --- main ---
! --comment>>
!git でバージョン管理
!Two-wayビームパターンの計算プログラム

program main
    use parameter_antenna
    use omp_lib
    implicit none
    integer :: i, j
    integer :: length, status, idx, local_scatnum
    real(8) :: t1, t2
    character(50) :: log_file, log_dir
    character(:), allocatable :: arg
    intrinsic :: command_argument_count, get_command_argument

    !
    ! コマンドライン引数を全て表示する．
    !
    do i = 0, command_argument_count()
        call get_command_argument(i, length = length, status = status)
        if (status == 0) then
        !
        ! 引数の文字長で "arg" をアロケートする．
        !
        allocate (character(length) :: arg)
        call get_command_argument(i, arg, status = status)
        if (status == 0) then
            !
            ! 取得した引数を表示する．
            !
            if (i == 0) then
                print *, 'Command = "', arg, '"'
            else
                print *, 'Argument', '= "', arg, '"'
                !
                !get_calctime
                !
                read(arg,*) idx
            end if
        end if
        deallocate (arg)
        end if
        !
        ! エラーの場合はエラーメッセージを表示する．
        !
        if (status /= 0) print *, 'Error', status, 'on argument', i
    end do

    !
    !call TX antenna_position
    !
!#if SUBARRAY_MODE==1
!    call antenna_position_allgrp
!#elif SUBARRAY_MODE==2
!    call antenna_position_sa
!#else
!    print*, "エラー：サブアレイの観測モードが適切に設定されていません。"
!    stop
!#endif
    call antenna_position_allgrp

    !
    !call RX antenna_position
    !
    call rx_antenna_position(100d0)

    !
    !散乱シミレーション
    !
    local_scatnum = int(num_scat/process)
    t1 = omp_get_wtime()
    call scattering_simulation_devided(idx, local_scatnum)
    t2 = omp_get_wtime()

    !
    !計算時間測定
    !
    write(log_dir, '("log_3Dscat_time/")')
    call makedirs(log_dir)
    write(log_file, '("idx", i3.3, "log.dat")') idx
    log_file = trim(log_dir)//trim(log_file)
    open(idx, file=log_file, status='replace')
    write(idx,*) t2-t1
    close(idx)

    !
    !antenna position
    !
    !print*, "==== antenna chech ===="
    !print*, "number of active radar modules", ACTIVE_GRPNUM
    !print*, "number of active radar antenna elements", ACTIVE_INGRP
    open(10, file='tx_antpos.dat', status='replace')
    do i = 1, ACTIVE_GRPNUM
       do j = 1, ANT_INGRP
          write(10,"(2(f6.1),i3)") ANT_XPOS(i,j), ANT_YPOS(i,j), i
       end do
    end do
    close(10)
    open(10, file='rx_antpos.dat', status='replace')
    do i = 1, RX_ANTNUM
        write(10,"(2f6.1,i3)") RX_XPOS(i), RX_YPOS(i), i
    end do
    close(10)

end program main
