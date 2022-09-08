#!/bin/bash

###プリプロセッサの設定
#PREPROCESSOR_SUBARRAY: {1:(subarray=25), 2:(subarray=sa mode), 3:(antenna test mode)}
PREPROCESSOR_SUBARRAY=1

###モデルパラメータの設定
#MAINLOBE_ZENITH:       メインローブの天頂角[度]
#MAINLOBE_AZIMUTH:      メインローブの方位角[度]
#THETA_X, THETA_Y:      観測体積の天頂角の大きさ[度]
#VX, VY:                観測体積内の水平面内の格子分割数
MAINLOBE_ZENITH="0d0"
MAINLOBE_AZIMUTH="0d0"
THETA_X="45d0"
THETA_Y="45d0"
VX=10000
VY=1000

###散乱シミュレーションのパラメータ設定
#DEG: モデルの回転角(度)
#RNG: 観測レンジ(m)
#ZCL: 計算領域の厚み(m)
#HVL: 水平平均風速(m/s)
#WVL: 上昇流(m/s)
#VSG: 速度分散(m/s) ---等方性乱流を仮定しているので速度分散は座標回転によらない
#NSC: 散乱体数(個)
#MTS: サンプリング回数
DEG="0d0"
RNG="1000d0"
ZCL="500d0"
HVL="0d0"
WVL="5d0"
VSG="1d0"
NSC="1e8"
MTS="2**12"

###散乱シミュレーションの計算パラメータ設定
#PROC:              ノード並列をする際に立ち上げるプロセス数
#OMP_NUM_THREADS:   OPENMPのスレッド数(プロセスあたりのコア数)
PROC="10"
OMP_NUM_THREADS="68"



###shell script内部の変数
ORIGINAL_PARAMETER_FILE="./parameter_gen.f90"
PRODUCED_PARAMETER_FILE="./parameter.f90"

###モデルパラメータの整合性の確認
#観測モードとアクティブサブアレイ数の確認
#ACTIVE_GRPNUM: 送信を行うサブアレイの個数
if test ${PREPROCESSOR_SUBARRAY} -eq 1; then
    ACTIVE_GRPNUM=25
    ACTIVE_INGRP=19
elif test ${PREPROCESSOR_SUBARRAY} -eq 2; then
    ACTIVE_GRPNUM=3
    ACTIVE_INGRP=19
elif test ${PREPROCESSOR_SUBARRAY} -eq 3; then
    ACTIVE_GRPNUM=2
    ACTIVE_INGRP=1
else
    cat << END
    [エラー!] サブアレイのプリプロセッサ設定はPREPROCESSOR_SUBARRAY={1:(subarray=25), 2:(subarray=sa mode)}のいずれかです。
    現在、PREPROCESSOR_SUBARRAY=${PREPROCESSOR_SUBARRAY}に設定されています。
END
    exit
fi
#格子数の確認
if [ `expr ${VX} % 2` == 1 ]; then
    echo "[エラー!]integer parameter: VX=${VX}, 偶数を設定してください。"
    exit
fi
if [ `expr ${VY} % 2` == 1 ]; then
    echo "[エラー!]integer parameter: VY=${VY}, 偶数を設定してください。"
    exit
fi

###ログディレクトリ作成
if [ ! -d "log" ]; then
    mkdir "log"
fi

###パラメータファイルの生成
if [ ! -e ${ORIGINAL_PARAMETER_FILE} ]; then
    echo "[エラー!]program file: ${ORIGINAL_PARAMETER_FILE} が存在しません"
    exit
else
    cp ${ORIGINAL_PARAMETER_FILE} ${PRODUCED_PARAMETER_FILE}
    echo "${PRODUCED_PARAMETER_FILE} is generated!"
    sed -i -e "s/__VX__/${VX}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__VY__/${VY}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__THETA_X__/${THETA_X}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__THETA_Y__/${THETA_Y}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__MAINLOBE_ZENITH__/${MAINLOBE_ZENITH}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__MAINLOBE_AZIMUTH__/${MAINLOBE_AZIMUTH}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__ACTIVE_GRPNUM__/${ACTIVE_GRPNUM}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__ACTIVE_INGRP__/${ACTIVE_INGRP}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__deg__/${DEG}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__rrange__/${RNG}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__zscale__/${ZCL}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__hvel_avr__/${HVL}/"  ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__w_avr__/${WVL}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__vel_sigma__/${VSG}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__num_scat__/${NSC}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__t_stepmax__/${MTS}/" ${PRODUCED_PARAMETER_FILE}
    sed -i -e "s/__process__/${PROC}/" ${PRODUCED_PARAMETER_FILE}
fi

###コンパイル
if ftn      -DSUBARRAY_MODE=${PREPROCESSOR_SUBARRAY} \
            -O3 \
            ${PRODUCED_PARAMETER_FILE} \
            subroutine.f90 \
            main_scat.f90
then
    rm *.mod
    #rm ${PRODUCED_PARAMETER_FILE}
    cat << END
    ################################################

    [コンパイル結果]
    COMPILE PROGRAM >>>>> SUCCEEDED!

    [シェルスクリプトで設定されたパラメータの表示]
    ##OPENMP THREADS
    OMP_NUM_THREADS=${OMP_NUM_THREADS}

    ##PREPROCESSOR
    -------------------------------------------------------
    -DBEAM_WRITE=${PREPROCESSOR_BEAM}           |{1:(ビームパターンの出力), other int:(出力なし)}
    -DSUBARRAY_MODE=${PREPROCESSOR_SUBARRAY}    |{1:(subarray=25), 2:(subarray=sa mode)}
    -DKOUSHI_MODE=${PREPROCESSOR_KOUSHI}        |{0:(等間隔格子),　1:(不等間隔格子)}
    -DXCF_MODE=${PREPROCESSOR_XCF}              |{0:(サブアレイごとのフィールド相関関数),　1:(全サブアレイ合成), 2:(サブアレイを合成)　2:(0,1,2の両方を計算)}
    -------------------------------------------------------

    ##${PRODUCED_PARAMETER_FILE}の設定
    -------------------------------------------------------
    parameter                   |value
    -------------------------------------------------------
    X軸方向格子数:              |${VX}
    Y軸方向格子数:              |${VY}
    X軸ビーム角:                |${THETA_X} [度]
    Y軸ビーム角:                |${THETA_Y} [度]
    アクティブなサブアレイ数:   |${ACTIVE_GRPNUM}
    アクティブなアンテナ素子:   |${ACTIVE_INGRP}
    メインローブ天頂角:         |${MAINLOBE_ZENITH} [度]
    メインローブ方位角:         |${MAINLOBE_AZIMUTH} [度]
    -------------------------------------------------------
    モデルの回転方向:           |${DEG} [度]
    観測レンジ:                 |${RNG} [m]
    計算領域の厚み:             |${ZCL} [m]
echo "hello $a"
    水平平均風速:               |${HVL} [m/s]
    上昇流:                     |${WVL} [m/s]
    速度分散:                   |${VSG} [m/s]
    散乱ターゲット数:           |${NSC} [個]
    サンプリング回数:           |${MTS} [回]
    プロセス数:                 |${PROC} [個]
    散乱ターゲット数/proc:      |${NSC}/${PROC} [個]
    -------------------------------------------------------
    OpenMP スレッド数:          |${OMP_NUM_THREADS}
    ################################################
END
fi

#multiprocess generation
OUTFILE="multiprocess.sh"
{
echo "#!/bin/bash"
echo "case \$ALPS_APP_PE in"
} > ${OUTFILE}
for ((i=0; i<${PROC}; i++))
do
echo "$i) ./a.out ${i} ;;" >> ${OUTFILE}
done
echo "esac" >> ${OUTFILE}
cat ${OUTFILE}

#execution script
OUTFILE="executoin_sysA.sh"
cat > ${OUTFILE} << END
#!/bin/bash
#que の指定は qstat -q で投入可能なのもを調べる
#============ PBS Options ============
#QSUB -q gr20001a
#QSUB -ug gr20001
#QSUB -W 60:00
#QSUB -A p=${PROC}:t=${OMP_NUM_THREADS}:c=${OMP_NUM_THREADS}:m=3G
#QSUB -M tamura.ryosuke.65w@st.kyoto-u.ac.jp
#QSUB -m b
#QSUB -m e
#============ Shell Script ============

aprun -n \${QSUB_PROCS} -d \${QSUB_THREADS} -N \${QSUB_PPN} -cc none -b sh ./multiprocess.sh
END
cat ${OUTFILE}
