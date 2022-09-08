! -------fortran program------- !
! --- subroutine ---
! --comment>>
!Two-wayビームパターンの計算プログラムのサブルーチン

!
!検査体積の格子座標生成を行うサブルーチン
!
subroutine observation_volume(range_local)
    use parameter_antenna
    implicit none
    real(8), intent(in) :: range_local
    real(8) :: xscale, yscale, x_center, y_center
    integer :: i
    xscale = 2d0 * range_local * tan(THETA_X/180d0*pi)
    yscale = 2d0 * range_local * tan(THETA_Y/180d0*pi)
    x_center = range_local * sin(mainlobe_zenith_deg/180d0*pi) * cos(mainlobe_azimuth_deg/180d0*pi)
    y_center = range_local * sin(mainlobe_zenith_deg/180d0*pi) * sin(mainlobe_azimuth_deg/180d0*pi)
    do i = 1, xgrid_num
        VH_X(i) = xscale/(xgrid_num-1)*(i-1) - xscale/2d0 + x_center
    enddo
    do i = 1, ygrid_num
        VH_Y(i) = yscale/(ygrid_num-1)*(i-1) - yscale/2d0 + y_center
    enddo
    !
    !check
    !
    print*, '========== observation volume grid generated ========='
    print"(a20,f10.2,a8)", 'radar range', range_local, '[m]'
    print"(a20,f10.2,a8)", 'main lobe zenith', mainlobe_zenith_deg, '[deg]'
    print"(a20,f10.2,a8)", 'main lobe azimuth', mainlobe_azimuth_deg, '[deg]'
    print"(a20,f10.2,a8)", 'center of volume x', x_center, '[m]'
    print"(a20,f10.2,a8)", 'center of volume y', y_center, '[m]'
    print*, '======================================================'
end subroutine observation_volume

!
!キャリブレーション
!
subroutine calibration
    use parameter_antenna
    implicit none
    real(8) avr_x, avr_y, avr_z
    avr_x = sum(ANT_XPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
    avr_y = sum(ANT_YPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
    avr_z = sum(ANT_ZPOS(1:ACTIVE_GRPNUM,:))/ALL_ANTNUM
    open(9000, file=logfile, position="append")
    write(9000,*) '================ check  calibration ========================'
    write(9000,*) 'Antenna(Totall)', ALL_ANTNUM
    write(9000,*) 'Center of garvity of antenna is (x,y,z)::'
    write(9000,"(3(f10.2,'m'))") avr_x, avr_y, avr_z
    close(9000)
end subroutine calibration


!
!2way beam pattern を計算する関数
!

!local_range:   レーダーのレンジ
!x_targ:        観測体積X座標
!y_targ:        観測体積Y座標
!z_targ:        観測体積Z座標
function func_reflex(local_range, x_targ, y_targ, z_targ) result(reflex)
    use parameter_antenna
    implicit none
    !
    !function
    !
    interface
        !
        !ガウス関数
        !
        function gaussian(x, mean, var) result(output)
            real(8), intent(in) :: x, mean, var
            real(8) :: output
        end function gaussian
    end interface
    real(8), intent(in) :: local_range(RX_ANTNUM), x_targ, y_targ, z_targ
    !
    !local 変数
    !
    complex(8) :: reflex(RX_ANTNUM)
    complex(8), parameter :: img = (0d0, 1d0)
    integer :: sub_tx
    integer :: elm_tx, elm_rx
    real(8) :: tx_targ, targ_rx
    real(8) :: phase
    !
    !initialize
    !
    reflex = (0d0, 0d0)
    !
    !calculation of signal
    !
    do elm_rx = 1, RX_ANTNUM
        targ_rx = sqrt( (x_targ - RX_XPOS(elm_rx))**2&
                       +(y_targ - RX_YPOS(elm_rx))**2&
                       +(z_targ - RX_ZPOS(elm_rx))**2 )
        do sub_tx = 1, ACTIVE_GRPNUM
            do elm_tx = 1, ACTIVE_INGRP
                tx_targ = sqrt( (x_targ - ANT_XPOS(sub_tx, elm_tx))**2&
                               +(y_targ - ANT_YPOS(sub_tx, elm_tx))**2&
                               +(z_targ - ANT_ZPOS(sub_tx, elm_tx))**2 )
                phase = kk*(tx_targ+targ_rx-ANT_PHASE(sub_tx, elm_tx))
                phase = mod(phase, 2d0*pi)
                reflex(elm_rx) = reflex(elm_rx) &
                               + exp( -img*phase )*rrange**2/(targ_rx*tx_targ) &
                               * gaussian(tx_targ+targ_rx, local_range(elm_rx), 75d0/sqrt(2d0*log(2d0)))
            end do
        end do
    end do
    return
end function func_reflex



!
!散乱シミレーションの実行
!

!idx:               計算プロセス番号
!local_scatnum:     ノード内で計算するターゲット数
subroutine scattering_simulation_devided(idx,local_scatnum)
    use parameter_antenna
    !$ use omp_lib
    implicit none
    integer, intent(in) :: idx, local_scatnum
    !
    !function
    !
    interface
        function func_reflex(local_range, x_targ, y_targ, z_targ) result(reflex)
            use parameter_antenna
            real(8), intent(in) :: local_range(RX_ANTNUM), x_targ, y_targ, z_targ
            complex(8) :: reflex(RX_ANTNUM)
        end function func_reflex
        function gaussian(x, mean, var) result(output)
            real(8), intent(in) :: x, mean, var
            real(8) :: output
        end function gaussian
    end interface
    real(8) :: xscale, yscale
    real(8) :: x_center, y_center, z_center
    real(8) :: mainlobe_zenith_rad, mainlobe_azimuth_rad
    real(8) :: local_range(RX_ANTNUM)
    !
    ! --- causion: local_scatnum must be integer --- !
    !
    real(8), allocatable :: XS(:), YS(:), ZS(:)
    real(8), allocatable :: VELX(:), VELY(:), VELZ(:)
    complex(8), allocatable :: SIGNAL_TMP(:,:)
    complex(8), parameter :: img = (0, 1.0)     !虚数単位、散乱強度
    integer :: t_step, seedsize, c, k, rx_idx
    integer, allocatable :: seed(:)
    real(8) :: vx_avr, vy_avr, deg_rad
    character(500) :: dirname, filename
    !
    !受信アンテナの重み付け距離
    !
    mainlobe_zenith_rad = mainlobe_zenith_deg/180d0*pi
    mainlobe_azimuth_rad = mainlobe_azimuth_deg/180d0*pi
    x_center = rrange*sin(mainlobe_zenith_rad)*cos(mainlobe_azimuth_deg)
    y_center = rrange*sin(mainlobe_zenith_rad)*sin(mainlobe_azimuth_deg)
    z_center = rrange*cos(mainlobe_zenith_rad)
    do rx_idx = 1, RX_ANTNUM
        ! local_range(rx_idx) = sqrt( (RX_XPOS(rx_idx)-x_center)**2 &
        !                           + (RX_YPOS(rx_idx)-y_center)**2 &
        !                           + (RX_ZPOS(rx_idx)-z_center)**2 &
        !                       ) + rrange
        ! print *, 'rx idx', rx_idx, 'range', local_range(rx_idx)
        local_range(rx_idx) = 2d0 * rrange
    end do
    !
    !region
    !
    xscale = rrange * tan(THETA_X/180.0*pi)
    yscale = rrange * tan(THETA_Y/180.0*pi)
    !
    !allocation
    !
    allocate(XS(local_scatnum))
    allocate(YS(local_scatnum))
    allocate(ZS(local_scatnum))
    allocate(VELX(local_scatnum))
    allocate(VELY(local_scatnum))
    allocate(VELZ(local_scatnum))
    allocate(SIGNAL_TMP(RX_ANTNUM,local_scatnum))
    !
    !randam setting
    !
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    call random_seed(get=seed(:))
    call system_clock(count=c)
    seed(1) = c + idx
    ! seed = 0
    call random_seed(put=seed(:))
    !
    !initial position of scatters
    !
    call random_number(XS)
    call random_number(YS)
    call random_number(ZS)
    XS = (2d0*XS-1)*xscale
    YS = (2d0*YS-1)*yscale
    ZS = (2d0*ZS-1)*zscale + rrange
    !
    !log
    !
    write(dirname,'("signal/mainlobe_zenith", i2.2, "mainlobe_azimuth", i3.3, "_TXsubarray", i2.2, "_RXelements", i2.2, "/")' &
                    &) int(mainlobe_zenith_deg), int(mainlobe_azimuth_deg), ACTIVE_GRPNUM, RX_ANTNUM
    call makedirs(dirname)
    write(filename, '("particle", e7.1, "dlength", i5.5, "range", i5.5, &
        & "wdirec", i3.3, "hvel", i3.3, f0.2, "wvel", i2.2, f0.2, "wsigma", i2.2, f0.2, "_idx", i6.6, ".dat")') &
        & real(num_scat), t_stepmax, int(rrange), int(windvector_deg), &
        & int(vel_avr), vel_avr-int(vel_avr), int(w_avr), w_avr-int(w_avr), int(vel_sigma), vel_sigma-int(vel_sigma), idx
    filename = trim(dirname)//trim(filename)
    !
    !mean wind velocity
    !
    deg_rad = windvector_deg/180d0*pi
    vx_avr = vel_avr*cos(deg_rad)
    vy_avr = vel_avr*sin(deg_rad)
    !
    !wind velocity of each scatters
    !
    call box_muller(local_scatnum, VELX, vx_avr, vel_sigma)
    call box_muller(local_scatnum, VELY, vy_avr, vel_sigma)
    call box_muller(local_scatnum, VELZ, w_avr, vel_sigma)
    !
    !radar signal
    !
    open(idx+500, file=filename, status='replace')
    do t_step = 1, t_stepmax
        !
        !posintion of scatters
        !
        XS = spltime * VELX + XS
        YS = spltime * VELY + YS
        ZS = spltime * VELZ + ZS
        !
        !boudary conditions
        !
        XS = modulo(XS + xscale, 2*xscale) - xscale
        YS = modulo(YS + yscale, 2*yscale) - yscale
        ZS = modulo(ZS + (zscale - rrange), 2*zscale) - (zscale - rrange)

        SIGNAL_TMP = (0d0, 0d0)
        !
        !散乱の計算
        !
        !$OMP parallel
        !$OMP do
        do k = 1, local_scatnum
            SIGNAL_TMP(:,k) = func_reflex(local_range, XS(k), YS(k), ZS(k)) + SIGNAL_TMP(:,k)
        end do
        !$OMP enddo
        !$OMP end parallel
        !
        !check
        !
        write(idx+500,*) (t_step-1)*spltime, real(sum(SIGNAL_TMP,dim=2)), aimag(sum(SIGNAL_TMP,dim=2))
        ! write(*,*) 'time step', t_step
    end do
    close(idx+500)
    !
    !deallocation
    !
    deallocate(XS, YS, ZS)
    deallocate(VELX, VELY, VELZ)
    deallocate(seed)
    deallocate(SIGNAL_TMP)
end subroutine scattering_simulation_devided


!
!Two-wayパターンの計算
!
subroutine two_way_pattern_openmp
    use parameter_antenna
    !$ use omp_lib
    implicit none
    !
    !functoin for calculation
    !
    interface
        !
        !gaussian function
        !
        function gaussian(x, mean, var) result(output)
            real(8), intent(in) :: x, mean, var
            real(8) :: output
        end function gaussian
        !
        !beam pattern
        !
        function func_reflex(local_range, x_targ, y_targ, z_targ) result(reflex)
            use parameter_antenna
            real(8), intent(in) :: local_range(RX_ANTNUM), x_targ, y_targ, z_targ
            complex(8) :: reflex(RX_ANTNUM)
        end function func_reflex
    end interface
    integer:: xgrid_idx, ygrid_idx, rx_idx
    real(8) :: x_center, y_center, z_center
    real(8) :: mainlobe_zenith_rad, mainlobe_azimuth_rad, local_range(RX_ANTNUM)
    complex(8), allocatable :: E_tmp(:,:,:)
    complex(8), parameter :: img = (0d0, 1d0)     !虚数単位、散乱強度

    !two_way_pattern記録
    character(500):: binary_file, header_file, dirname
    write(dirname,'("beampattern/mainlobe_zenith", i2.2, "mainlobe_azimuth", i3.3, "_subarray", i2.2,&
                    &"_zenithX",  i3.3, "VX", i7.7, &
                    &"_zenithY",  i3.3, "VY", i7.7, "/")' &
                    &) int(mainlobe_zenith_deg), int(mainlobe_azimuth_deg), ACTIVE_GRPNUM, &
                    int(THETA_X), xgrid_num, int(THETA_Y), ygrid_num
    call makedirs(dirname)

    allocate(E_tmp(RX_ANTNUM, xgrid_num, ygrid_num))
    E_tmp(:,:,:) = 0d0
    call observation_volume(rrange)
    mainlobe_zenith_rad = mainlobe_zenith_deg/180d0*pi
    mainlobe_azimuth_rad = mainlobe_azimuth_deg/180d0*pi
    x_center = rrange*sin(mainlobe_zenith_rad)*cos(mainlobe_azimuth_deg)
    y_center = rrange*sin(mainlobe_zenith_rad)*sin(mainlobe_azimuth_deg)
    z_center = rrange*cos(mainlobe_zenith_rad)
    !
    !受信アンテナの重み付け距離
    !
    do rx_idx = 1, RX_ANTNUM
        ! local_range(rx_idx) = sqrt( (RX_XPOS(rx_idx)-x_center)**2 &
        !                           + (RX_YPOS(rx_idx)-y_center)**2 &
        !                           + (RX_ZPOS(rx_idx)-z_center)**2 &
        !                       ) + rrange
        ! print *, 'rx idx', rx_idx, 'range', local_range(rx_idx)
        local_range(rx_idx) = 2d0 * rrange
    end do
    !$OMP parallel
    !$OMP do
    do ygrid_idx = 1, ygrid_num
        do xgrid_idx = 1, xgrid_num
            !
            !two way beam pattern
            !
            E_tmp(:,xgrid_idx,ygrid_idx) = func_reflex(local_range, VH_X(xgrid_idx), VH_Y(ygrid_idx), z_center)
        end do
    end do
    !$OMP enddo
    !$OMP end parallel
    print*, '========== beam pattern calculation ========='
    print"(a20,f10.2,a8)", 'radar range', rrange, '[m]'
    print*, '=============================================='

    !
    !beam pattern data
    !
    write(binary_file,'("binary_range", i6.6, ".bin")') &
                   & int(rrange)
    binary_file = trim(dirname)//trim(binary_file)
    open(8000, file=binary_file, form='unformatted', status='replace')
    write(8000) real(E_tmp), aimag(E_tmp)
    close(8000)
    !
    !data header
    !
    write(header_file,'("header_range", i6.6, ".txt")')&
                   & int(rrange)
    header_file = trim(dirname)//trim(header_file)
    open(7000, file=header_file, status='replace')
    write(7000, '(a20,a)') 'bfile ', trim(binary_file)
    write(7000, '(a20,i8)') 'subarrynum_tx', ACTIVE_GRPNUM
    write(7000, '(a20,i8)') 'subarrynum_rx', RX_ANTNUM
    write(7000, '(a20,i8)') 'x_gridnum', xgrid_num
    write(7000, '(a20,i8)') 'y_gridnum', ygrid_num
    write(7000, '(a20,f8.1)') 'x_min', VH_X(1)
    write(7000, '(a20,f8.1)') 'x_max', VH_X(xgrid_num)
    write(7000, '(a20,f8.1)') 'y_min', VH_Y(1)
    write(7000, '(a20,f8.1)') 'y_max', VH_Y(ygrid_num)
    write(7000, '(a20,f8.1)') 'radar_range', rrange
    write(7000, '(a20,f8.1)') 'mainlobe_zenith', mainlobe_zenith_deg
    write(7000, '(a20,f8.1)') 'mainlobe_azimuth', mainlobe_azimuth_deg
    close(7000)
    deallocate(E_tmp)
end subroutine two_way_pattern_openmp


!
!ガウス関数
!
function gaussian(x, mean, var) result(output)
    implicit none
    real(8), intent(in) :: x, mean, var
    real(8) :: output
    output = exp( -(x-mean)**2/(2.0*var**2) )
    return
end function gaussian


!
!boxmuller method
!正規乱数の生成
!
subroutine box_muller(nsample, rnd, mean, std)
    implicit none
    integer, intent(in) :: nsample
    real(8), intent(in) :: mean, std
    real(8), intent(out) :: rnd(nsample)
    real(8), allocatable :: x(:), y(:), r(:), t(:)
    real(8), parameter :: pi = 3.141592653589793
    integer :: i, seedsize, c
    integer, allocatable :: seed(:)

    !
    !allocation
    !
    allocate(x(nsample), y(nsample), r(nsample), t(nsample))

    !
    !random setting
    !
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    call random_seed(get=seed(:))

    !
    !random depends on time
    !
    call system_clock(count=c)
    seed(1) = c
    ! seed = 1
    call random_seed(put=seed(:))

    !
    !randam variable
    !
    call random_number(x)
    call random_number(y)

    !
    !normal randam variable N(0,1)
    !
    r = sqrt(-2d0 * log(x))
    t = 2d0 * pi * y
    rnd = r * cos(t)
    !rnd = r * sin(t)
    !requred randam variable N(mean, std)
    do i = 1, nsample
        rnd(i) = rnd(i)*std + mean
    end do

    !
    !deallocation
    !
    deallocate(seed)
end subroutine box_muller


!
!directory を作成するサブルーチン
!
subroutine makedirs(outdir)
    character(len=*), intent(in) :: outdir
    character(len=900) command
    write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
    !write(*, *) trim(command)
    call system(command)
end subroutine makedirs
