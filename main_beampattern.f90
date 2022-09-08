! -------fortran program------- !
! --- main ---
! --comment>>
!git でバージョン管理
!Two-wayビームパターンの計算プログラム

program main
    use parameter_antenna
    implicit none
    integer :: i, j

    !
    !call TX antenna_position
    !
#if SUBARRAY_MODE==1
    call antenna_position_allgrp
#elif SUBARRAY_MODE==2
    call antenna_position_sa
#else
    print*, "エラー：サブアレイの観測モードが適切に設定されていません。"
    stop
#endif

    !
    !call RX antenna_position
    !
    call rx_antenna_position(100d0)

    !
    !ビームパターン
    !
    call two_way_pattern_openmp

    !
    !antenna position
    !
    open(10, file='tx_antpos.dat', status='replace')
    do i = 1, ACTIVE_GRPNUM
       do j = 1, ANT_INGRP
          write(10,"(2(f7.1),i3)") ANT_XPOS(i,j), ANT_YPOS(i,j), i
       end do
    end do
    close(10)
    open(10, file='rx_antpos.dat', status='replace')
    do i = 1, RX_ANTNUM
        write(10,"(2f7.1,i3)") RX_XPOS(i), RX_YPOS(i), i
    end do
    close(10)
end program main
