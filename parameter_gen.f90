! -------fortran program------- !
! --- parameter ---
! --comment>>

!
!パラメータの設定
!
module parameter_antenna
    implicit none
    !
    !物理パラメータ
    !
    !cc:            光速 [m/s]
    !pi:            円周率
    !ff:            MUレーダーの中心周波数 [Hz]
    !kk:            MUレーダーの波数 [m^(-1)]
    real(8), parameter :: cc = 299792458d0
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: ff = 46.5*1e6
    real(8), parameter :: kk = 2d0*pi*ff/cc

    !
    !観測対象体積のパラメータ
    !
    !VX, VY:                観測体積内の水平面内の格子分割数(モジュール内で指定)
    !VH_X, VH_Y:            観測体積の水平格子座標配列
    !VLX, VLY:              観測体積の水平面の長さ [m]
    !THETA_X, THETA_Y:      X,Y方向のレーダービーム天頂角 [度]
    !mainlobe_zenith_deg:   メインローブの天頂角 [度]
    !mainlobe_azimuth_deg:  メインローブの方位角 [度]
    !RADAR_RANGE:           観測レンジ [m]
    !VEL_DEG:               モデルの回転角度(半時計回りを正)[度]
    !TILT:                  水平風と上昇流の成す角度(上向きを正)[度]
    integer, parameter :: xgrid_num = __VX__, ygrid_num = __VY__
    real(8) :: VH_X(xgrid_num+1), VH_Y(ygrid_num+1)
    real(8) :: VLX, VLY
    real(8), parameter :: THETA_X = __THETA_X__, THETA_Y = __THETA_Y__
    real(8), parameter :: mainlobe_zenith_deg = __MAINLOBE_ZENITH__
    real(8), parameter :: mainlobe_azimuth_deg = __MAINLOBE_AZIMUTH__

    !
    !散乱シミュレーションのパラメター
    !
    !windvector_deg:    散乱体の進む方向
    !rrange:            観測レンジ
    !vel_avr:           水平平均速度
    !w_avr:             上昇流
    !vel_sigma:         速度分散 ---等方性乱流を仮定しているので速度分散は座標回転によらない
    real(8), parameter :: windvector_deg = __deg__
    real(8), parameter :: rrange = __rrange__
    real(8), parameter :: vel_avr = __hvel_avr__
    real(8), parameter :: w_avr = __w_avr__
    real(8), parameter :: vel_sigma = __vel_sigma__
    integer, parameter :: num_scat = __num_scat__
    integer, parameter :: t_stepmax = __t_stepmax__
    integer, parameter :: process = __process__
    real(8), parameter :: spltime = 0.032d0
    real(8), parameter :: zscale = __zscale__

    !
    !アンテナの設定パラメータ
    !
    !ANT_INGRP:     1つのサブアレイを構成するアンテナ本数
    !ALL_GRPNUM:    サブアレイの総数
    !ACTIVE_GRPNUM: 送信を行うサブアレイの個数
    !ACTIVE_INGRP:  送信を行うサブアレイ中のアンテナ本数
    !ALL_ANTNUM:    送信を行うアンテナの本数
    !BEAM_MODE:     ビームフォーミング方向, {0:(0,0) 1:N(0,10) 2:E(90,10), 3:S(180,10), 4:W(270,10)}
    !ANT_*POS:      アンテナの位置座標配列 [m]
    !ANT_PHASE:     ビームフォーミングに対応した加算位相配列 [m]
    !RX_ANTNUM:     受信アンテナの本数
    !RX_*POS:       受信アンテナの座標
    integer, parameter :: ANT_INGRP = 19, ALL_GRPNUM = 25
    integer, parameter :: ACTIVE_GRPNUM = __ACTIVE_GRPNUM__
    integer, parameter :: ACTIVE_INGRP = __ACTIVE_INGRP__
    integer :: ALL_ANTNUM = ANT_INGRP * ACTIVE_GRPNUM
    real(8), dimension(ALL_GRPNUM,ANT_INGRP) :: ANT_XPOS, ANT_YPOS, ANT_ZPOS
    real(8), dimension(ALL_GRPNUM,ANT_INGRP) :: ANT_PHASE
    integer, parameter :: RX_ANTNUM = 5
    real(8), dimension(RX_ANTNUM) :: RX_XPOS, RX_YPOS, RX_ZPOS

    !
    !その他の変数
    !
    !logfile:       計算のログを記録
    !kousi_flag:    不等間隔・等間隔を表すフラグ
    character(400):: logfile

    !
    !MUレーダアンテナの座標
    !
    real(8), parameter :: A1_X(19) = (/  8,  7,  6,  5,  5,  4,  4,  3,  3,  2,  2,  1,  0,  0, -1, -1, -2, -2, -3/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: A2_X(19) = (/  3,  3,  3,  2,  2,  2,  2,  1,  1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: A3_X(19) = (/  8,  8,  8,  7,  7,  7,  7,  6,  6,  6,  6,  6,  5,  5,  5,  5,  4,  4,  4/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: A4_X(19) = (/  5,  5,  5,  4,  4,  4,  4,  3,  3,  3,  3,  3,  2,  2,  2,  2,  1,  1,  1/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: B1_X(19) = (/ 13, 13, 13, 13, 12, 12, 12, 12, 11, 11, 11, 11, 11, 10, 10, 10,  9,  9,  9/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: B2_X(19) = (/ 10, 10, 10,  9,  9,  9,  9,  8,  8,  8,  8,  8,  7,  7,  7,  7,  6,  6,  6/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: B3_X(19) = (/ 12, 12, 12, 11, 11, 11, 11, 10, 10, 10, 10, 10,  9,  9,  9,  9,  8,  8,  8/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: B4_X(19) = (/  7,  7,  7,  6,  6,  6,  6,  5,  5,  5,  5,  5,  4,  4,  4,  4,  3,  3,  3/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: C1_X(19) = (/ 12, 12, 11, 11, 11, 10, 10, 10, 10,  9,  9,  8,  8,  8,  7,  7,  7,  6,  5/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: C2_X(19) = (/  9,  9,  9,  8,  8,  8,  8,  7,  7,  7,  7,  7,  6,  6,  6,  6,  5,  5,  5/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: C3_X(19) = (/  6,  6,  6,  5,  5,  5,  5,  4,  4,  4,  4,  4,  3,  3,  3,  3,  2,  2,  2/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: C4_X(19) = (/  4,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  1,  1,  1,  1,  0,  0,  0/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: D1_X(19) = (/  3,  2,  2,  1,  1,  0,  0, -1, -2, -2, -3, -3, -4, -4, -5, -5, -6, -7, -8/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: D2_X(19) = (/  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1, -1, -2, -2, -2, -2, -3, -3, -3/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: D3_X(19) = (/ -4, -4, -4, -5, -5, -5, -5, -6, -6, -6, -6, -6, -7, -7, -7, -7, -8, -8, -8/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: D4_X(19) = (/ -1, -1, -1, -2, -2, -2, -2, -3, -3, -3, -3, -3, -4, -4, -4, -4, -5, -5, -5/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: E1_X(19) = (/ -9, -9, -9,-10,-10,-10,-11,-11,-11,-11,-11,-12,-12,-12,-12,-13,-13,-13,-13/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: E2_X(19) = (/ -6, -6, -6, -7, -7, -7, -7, -8, -8, -8, -8, -8, -9, -9, -9, -9,-10,-10,-10/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: E3_X(19) = (/ -8, -8, -8, -9, -9, -9, -9,-10,-10,-10,-10,-10,-11,-11,-11,-11,-12,-12,-12/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: E4_X(19) = (/ -3, -3, -3, -4, -4, -4, -4, -5, -5, -5, -5, -5, -6, -6, -6, -6, -7, -7, -7/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: F1_X(19) = (/ -5, -6, -7, -7, -7,- 8, -8, -8, -9, -9,-10,-10,-10,-10,-11,-11,-11,-12,-12/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: F2_X(19) = (/ -5, -5, -5, -6, -6, -6, -6, -7, -7, -7, -7, -7, -8, -8, -8, -8, -9, -9, -9/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: F3_X(19) = (/ -2, -2, -2, -3, -3, -3, -3, -4, -4, -4, -4, -4, -5, -5, -5, -5, -6, -6, -6/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: F4_X(19) = (/  0,  0,  0, -1, -1, -1, -1, -2, -2, -2, -2, -2, -3, -3, -3, -3, -4, -4, -4/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: F5_X(19) = (/  2,  2,  2,  1,  1,  1,  1,  0,  0,  0,  0,  0, -1, -1, -1, -1, -2, -2, -2/) &
                                    * 4.5 * sqrt(3d0) / 2d0
    real(8), parameter :: A1_Y(19) = (/ 18, 19, 20, 21, 19, 20, 18, 21, 19, 22, 20, 21, 22, 20, 21, 19, 22, 20, 21/) * 4.5 / 2d0
    real(8), parameter :: A2_Y(19) = (/ 17, 15, 13, 18, 16, 14, 12, 19, 17, 15, 13, 11, 18, 16, 14, 12, 17, 15, 13/) * 4.5 / 2d0
    real(8), parameter :: A3_Y(19) = (/ 16, 14, 12, 17, 15, 13, 11, 18, 16, 14, 12, 10, 17, 15, 13, 11, 16, 14, 12/) * 4.5 / 2d0
    real(8), parameter :: A4_Y(19) = (/  9,  7,  5, 10,  8,  6,  4, 11,  9,  7,  5,  3, 10,  8,  6,  4,  9,  7,  5/) * 4.5 / 2d0
    real(8), parameter :: B1_Y(19) = (/  3,  1, -1, -3,  8,  6,  4,  2, 11,  9,  7,  5,  3, 14, 12, 10, 15, 13, 11/) * 4.5 / 2d0
    real(8), parameter :: B2_Y(19) = (/  8,  6,  4,  9,  7,  5,  3, 10,  8,  6,  4,  2,  9,  7,  5,  3,  8,  6,  4/) * 4.5 / 2d0
    real(8), parameter :: B3_Y(19) = (/  0, -2, -4,  1, -1, -3, -5,  2,  0, -2, -4, -6,  1, -1, -3, -5,  0, -2, -4/) * 4.5 / 2d0
    real(8), parameter :: B4_Y(19) = (/  1, -1, -3,  2,  0, -2, -4,  3,  1, -1, -3, -5,  2,  0, -2, -4,  1, -1, -3/) * 4.5 / 2d0
    real(8), parameter :: C1_Y(19) = (/ -6, -8, -7, -9,-11, -8,-10,-12,-14,-13,-15,-14,-16,-18,-15,-17,-19,-20,-21/) * 4.5 / 2d0
    real(8), parameter :: C2_Y(19) = (/ -7, -9,-11, -6, -8,-10,-12, -5, -7, -9,-11,-13, -6, -8,-10,-12, -7, -9,-11/) * 4.5 / 2d0
    real(8), parameter :: C3_Y(19) = (/-14,-16,-18,-13,-15,-17,-19,-12,-14,-16,-18,-20,-13,-15,-17,-19,-14,-16,-18/) * 4.5 / 2d0
    real(8), parameter :: C4_Y(19) = (/ -6, -8,-10, -5, -7, -9,-11, -4, -6, -8,-10,-12, -5, -7, -9,-11, -6, -8,-10/) * 4.5 / 2d0
    real(8), parameter :: D1_Y(19) = (/-21,-20,-22,-19,-21,-20,-22,-21,-20,-22,-19,-21,-18,-20,-19,-21,-20,-19,-18/) * 4.5 / 2d0
    real(8), parameter :: D2_Y(19) = (/-13,-15,-17,-12,-14,-16,-18,-11,-13,-15,-17,-19,-12,-14,-16,-18,-13,-15,-17/) * 4.5 / 2d0
    real(8), parameter :: D3_Y(19) = (/-12,-14,-16,-11,-13,-15,-17,-10,-12,-14,-16,-18,-11,-13,-15,-17,-12,-14,-16/) * 4.5 / 2d0
    real(8), parameter :: D4_Y(19) = (/ -5, -7, -9, -4, -6, -8,-10, -3, -5, -7, -9,-11, -4, -6, -8,-10, -5, -7, -9/) * 4.5 / 2d0
    real(8), parameter :: E1_Y(19) = (/-11,-13,-15,-10,-12,-14, -3, -5, -7, -9,-11, -2, -4, -6, -8,  3,  1, -1, -3/) * 4.5 / 2d0
    real(8), parameter :: E2_Y(19) = (/ -4, -6, -8, -3, -5, -7, -9, -2, -4, -6, -8,-10, -3, -5, -7, -9, -4, -6, -8/) * 4.5 / 2d0
    real(8), parameter :: E3_Y(19) = (/  4,  2,  0,  5,  3,  1, -1,  6,  4,  2,  0, -2,  5,  3,  1, -1,  4,  2,  0/) * 4.5 / 2d0
    real(8), parameter :: E4_Y(19) = (/  3,  1, -1,  4,  2,  0, -2,  5,  3,  1, -1, -3,  4,  2,  0, -2,  3,  1, -1/) * 4.5 / 2d0
    real(8), parameter :: F1_Y(19) = (/ 21, 20, 19, 17, 15, 18, 16, 14, 15, 13, 14, 12, 10,  8, 11,  9,  7,  8,  6/) * 4.5 / 2d0
    real(8), parameter :: F2_Y(19) = (/ 11,  9,  7, 12, 10,  8,  6, 13, 11,  9,  7,  5, 12, 10,  8,  6, 11,  9,  7/) * 4.5 / 2d0
    real(8), parameter :: F3_Y(19) = (/ 18, 16, 14, 19, 17, 15, 13, 20, 18, 16, 14, 12, 19, 17, 15, 13, 18, 16, 14/) * 4.5 / 2d0
    real(8), parameter :: F4_Y(19) = (/ 10,  8,  6, 11,  9,  7,  5, 12, 10,  8,  6,  4, 11,  9,  7,  5, 10,  8,  6/) * 4.5 / 2d0
    real(8), parameter :: F5_Y(19) = (/  2,  0, -2,  3,  1, -1, -3,  4,  2,  0, -2, -4,  3,  1, -1, -3,  2,  0, -2/) * 4.5 / 2d0

contains

    !
    !テスト用アンテナ座標
    !
    subroutine antenna_position_test
        implicit none
        real(8) avr_x, avr_y, avr_z
        ANT_XPOS = 0d0
        ANT_YPOS = 0d0
        ANT_ZPOS = 0d0
        ANT_XPOS(1,:) = (0);   ANT_YPOS(1,:) = (0)
        ANT_XPOS(2,:) = (3);   ANT_YPOS(2,:) = (0)
        !calibration
        avr_x = sum(ANT_XPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        avr_y = sum(ANT_YPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        avr_z = sum(ANT_ZPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        ANT_XPOS = ANT_XPOS - avr_x
        ANT_YPOS = ANT_YPOS - avr_y
        ANT_ZPOS = ANT_ZPOS - avr_z
        !
        !ビームフォーミングの遅延量(距離)
        !
        call phase_delay
    end subroutine antenna_position_test

    !
    !SA 解析用のデータ
    !
    subroutine antenna_position_sa
        implicit none
        real(8) avr_x, avr_y, avr_z

        ANT_XPOS = 0d0
        ANT_YPOS = 0d0
        ANT_ZPOS = 0d0
        ANT_XPOS(1,:) = F2_X;   ANT_YPOS(1,:) = F2_Y
        ANT_XPOS(2,:) = F3_X;   ANT_YPOS(2,:) = F3_Y
        ANT_XPOS(3,:) = F4_X;   ANT_YPOS(3,:) = F4_Y

        !calibration
        avr_x = sum(ANT_XPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        avr_y = sum(ANT_YPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        avr_z = sum(ANT_ZPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        ANT_XPOS = ANT_XPOS - avr_x
        ANT_YPOS = ANT_YPOS - avr_y
        ANT_ZPOS = ANT_ZPOS - avr_z

        !
        !ビームフォーミングの遅延量(距離)
        !
        call phase_delay

    end subroutine antenna_position_sa


    !アンテナ座標を設定するサブルーチン
    !全群解析の場合
    subroutine antenna_position_allgrp
        implicit none
        real(8) avr_x, avr_y, avr_z

        ANT_XPOS = 0d0
        ANT_YPOS = 0d0
        ANT_ZPOS = 0d0
        ANT_XPOS( 1,:) = A1_X;  ANT_YPOS( 1,:) = A1_Y
        ANT_XPOS( 2,:) = A2_X;  ANT_YPOS( 2,:) = A2_Y
        ANT_XPOS( 3,:) = A3_X;  ANT_YPOS( 3,:) = A3_Y
        ANT_XPOS( 4,:) = A4_X;  ANT_YPOS( 4,:) = A4_Y
        ANT_XPOS( 5,:) = B1_X;  ANT_YPOS( 5,:) = B1_Y
        ANT_XPOS( 6,:) = B2_X;  ANT_YPOS( 6,:) = B2_Y
        ANT_XPOS( 7,:) = B3_X;  ANT_YPOS( 7,:) = B3_Y
        ANT_XPOS( 8,:) = B4_X;  ANT_YPOS( 8,:) = B4_Y
        ANT_XPOS( 9,:) = C1_X;  ANT_YPOS( 9,:) = C1_Y
        ANT_XPOS(10,:) = C2_X;  ANT_YPOS(10,:) = C2_Y
        ANT_XPOS(11,:) = C3_X;  ANT_YPOS(11,:) = C3_Y
        ANT_XPOS(12,:) = C4_X;  ANT_YPOS(12,:) = C4_Y
        ANT_XPOS(13,:) = D1_X;  ANT_YPOS(13,:) = D1_Y
        ANT_XPOS(14,:) = D2_X;  ANT_YPOS(14,:) = D2_Y
        ANT_XPOS(15,:) = D3_X;  ANT_YPOS(15,:) = D3_Y
        ANT_XPOS(16,:) = D4_X;  ANT_YPOS(16,:) = D4_Y
        ANT_XPOS(17,:) = E1_X;  ANT_YPOS(17,:) = E1_Y
        ANT_XPOS(18,:) = E2_X;  ANT_YPOS(18,:) = E2_Y
        ANT_XPOS(19,:) = E3_X;  ANT_YPOS(19,:) = E3_Y
        ANT_XPOS(20,:) = E4_X;  ANT_YPOS(20,:) = E4_Y
        ANT_XPOS(21,:) = F1_X;  ANT_YPOS(21,:) = F1_Y
        ANT_XPOS(22,:) = F2_X;  ANT_YPOS(22,:) = F2_Y
        ANT_XPOS(23,:) = F3_X;  ANT_YPOS(23,:) = F3_Y
        ANT_XPOS(24,:) = F4_X;  ANT_YPOS(24,:) = F4_Y
        ANT_XPOS(25,:) = F5_X;  ANT_YPOS(25,:) = F5_Y

        !calibration
        avr_x = sum(ANT_XPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        avr_y = sum(ANT_YPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        avr_z = sum(ANT_ZPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
        ANT_XPOS = ANT_XPOS - avr_x
        ANT_YPOS = ANT_YPOS - avr_y
        ANT_ZPOS = ANT_ZPOS - avr_z

        !
        !ビームフォーミングの遅延量(距離)
        !
        call phase_delay

    end subroutine antenna_position_allgrp



    !
    !外付け受信アンテナの座標
    !
    subroutine rx_antenna_position(leng_d)
        implicit none
        real(8), intent(in) :: leng_d
        real(8) :: rad
        integer :: idx_rx
        RX_XPOS(1) = 0d0
        RX_YPOS(1) = 0d0
        RX_YPOS(1) = 0d0
        rad = 90d0/180d0 * pi
        do idx_rx = 2, RX_ANTNUM
            RX_XPOS(idx_rx) = leng_d * cos(rad*(idx_rx-2))
            RX_YPOS(idx_rx) = leng_d * sin(rad*(idx_rx-2))
            RX_ZPOS(idx_rx) = 0d0
        end do
    end subroutine


    !
    !ビームフォーミングの位相遅延量の計算
    !
    subroutine phase_delay
        implicit none
        integer :: sub, elm
        real(8) :: theta, phi
        theta = mainlobe_zenith_deg/180d0*pi
        phi = mainlobe_azimuth_deg/180d0*pi
        do sub = 1, ACTIVE_GRPNUM
            do elm = 1, ACTIVE_INGRP
                ANT_PHASE(sub, elm) = - sin(theta)*( ANT_XPOS(sub,elm)*cos(phi) + ANT_YPOS(sub,elm)*sin(phi) )
            end do
        end do
    end subroutine phase_delay


end module parameter_antenna
