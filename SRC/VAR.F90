!==================================================================================================
! 变量模块 VAR (Variable Module)
! 主要功能: 定义变量
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


MODULE VAR

  IMPLICIT NONE

  !------------------------------------------------------------------------------------------------
  !           定义子程序包控制变量
  !------------------------------------------------------------------------------------------------

  ! 控制 AMP 子程序包的启用
  LOGICAL :: L_AMP_BAS
  LOGICAL :: L_AMP_PET
  LOGICAL :: L_AMP_PCP

  ! 控制 CRP 子程序包的启用
  LOGICAL :: L_CRP_BAS
  LOGICAL :: L_CRP_PET
  LOGICAL :: L_CRP_INT

  ! 控制 IRR 子程序包的启用
  LOGICAL :: L_IRR_BAS
  LOGICAL :: L_IRR_FWA
  LOGICAL :: L_IRR_OWA

  ! 控制 SWF 子程序包的启用
  LOGICAL :: L_SWF_BAS
  LOGICAL :: L_SWF_STO
  LOGICAL :: L_SWF_CCF
  LOGICAL :: L_SWF_GBB
  LOGICAL :: L_SWF_GTB
  LOGICAL :: L_SWF_INF
  LOGICAL :: L_SWF_EVP
  LOGICAL :: L_SWF_RWU
  LOGICAL :: L_SWF_SSK
  LOGICAL :: L_SWF_TMA

  ! 控制 SHT 子程序包的启用
  LOGICAL :: L_SHT_BAS
  LOGICAL :: L_SHT_STO
  LOGICAL :: L_SHT_CCT
  LOGICAL :: L_SHT_GBB
  LOGICAL :: L_SHT_GTB
  LOGICAL :: L_SHT_TMA

  ! 控制 SST 子程序包的启用
  LOGICAL :: L_SST_BAS
  LOGICAL :: L_SST_STO
  LOGICAL :: L_SST_CCD
  LOGICAL :: L_SST_CCA
  LOGICAL :: L_SST_GBB
  LOGICAL :: L_SST_GTB
  LOGICAL :: L_SST_INF
  LOGICAL :: L_SST_TMA


  !------------------------------------------------------------------------------------------------
  !           定义输入文件变量
  !------------------------------------------------------------------------------------------------

  ! 基本数据输入文件
  INTEGER, PARAMETER :: U_IN_BAS = 100
  CHARACTER(LEN = 256) :: S_IN_BAS = 'BAS.IN'

  ! 土壤数据输入文件
  INTEGER, PARAMETER :: U_IN_SOL = 120
  CHARACTER(LEN = 256) :: S_IN_SOL

  ! 应力期数据输入文件
  INTEGER, PARAMETER :: U_IN_PER = 130
  CHARACTER(LEN = 256) :: S_IN_PER

  ! 气象数据输入文件
  INTEGER, PARAMETER :: U_IN_MET = 140
  CHARACTER(LEN = 256) :: S_IN_MET

  ! 作物数据输入文件
  INTEGER, PARAMETER :: U_IN_CRP = 150
  CHARACTER(LEN = 256) :: S_IN_CRP


  !------------------------------------------------------------------------------------------------
  !           定义日志文件和输出文件变量
  !------------------------------------------------------------------------------------------------

  ! 日志文件
  INTEGER, PARAMETER :: U_LOG = 200
  CHARACTER(LEN = 256) :: S_LOG = 'LOG.OUT'

  ! 土壤剖面数据输出文件
  INTEGER, PARAMETER :: U_OUT_LST = 300
  CHARACTER(LEN = 256) :: S_OUT_LST

  ! 应力期数据输出文件
  INTEGER, PARAMETER :: U_OUT_PER = 310
  CHARACTER(LEN = 256) :: S_OUT_PER

  ! 观测点数据输出文件
  INTEGER, PARAMETER :: U_OUT_OBS = 330
  CHARACTER(LEN = 256) :: S_OUT_OBS


  !------------------------------------------------------------------------------------------------
  !           定义常规设置变量
  !------------------------------------------------------------------------------------------------

  ! 程序名称
  CHARACTER(LEN = 80), DIMENSION(5), PARAMETER :: APPNAM = (/                          &
    '================================================================================', &
    '                              农田水热盐运移数值模型                              ', &
    '                                    FIRI-WHS                                    ', &
    '                            中国农业科学院农田灌溉研究所                          ', &
    '================================================================================'/)

  ! 项目名称, 限 100 个英文字符
  CHARACTER(LEN = 100) :: PRJNAM

  ! 运行模式
  ! IRUNMOD = 1, 正常模式, 出错时退出程序
  !         = 2, 安静模式, 出错时退出程序
  !         = 3, 安静模式, 出错时继续运行
  INTEGER :: IRUNMOD

  ! 错误代码的返回值
  ! IERRCOD = 0, 正常
  !         = 1, 土壤水分运动方程组系数错误
  !         = 2, 土壤水分运动方程组的解不收敛
  !         = 3, 土壤盐分运移方程组系数错误
  !         = 4, 土壤盐分运移方程组的解不收敛
  INTEGER :: IERRCOD

  ! 控制是否输入气象数据
  ! LMETEO = .TRUE., 是
  !        = .FALSE., 否
  LOGICAL :: LMETEO

  ! 控制是否输入作物数据
  ! LCROP = .TRUE., 是
  !       = .FALSE., 否
  LOGICAL :: LCROP

  ! 控制是否输入灌溉数据
  ! LIRRI = .TRUE., 是
  !       = .FALSE., 否
  LOGICAL :: LIRRI

  ! 控制是否模拟土壤热传导
  ! LHEAT = .TRUE., 是
  !       = .FALSE., 否
  LOGICAL :: LHEAT

  ! 控制是否模拟土壤盐分运移
  ! LSALT = .TRUE., 是
  !       = .FALSE., 否
  LOGICAL :: LSALT

  ! 控制是否模拟根系吸水
  ! LRWU = .TRUE., 是
  !      = .FALSE., 否
  LOGICAL :: LRWU

  ! 控制是否模拟源/汇项
  ! LSSK = .TRUE., 是
  !      = .FALSE., 否
  LOGICAL :: LSSK


  !------------------------------------------------------------------------------------------------
  !           定义模型规模变量
  !------------------------------------------------------------------------------------------------

  ! 土壤材料数量
  INTEGER :: NMAT

  ! 土壤单元体数量
  INTEGER :: NCEL

  ! 应力期数量
  INTEGER :: NPER

  ! 每个应力期的时长 (d)
  REAL, DIMENSION(:), ALLOCATABLE :: PERLEN

  ! 每个应力期的时段数量
  INTEGER, DIMENSION(:), ALLOCATABLE :: NSTP

  ! 每个应力期的时段步长增量因子
  REAL, DIMENSION(:), ALLOCATABLE :: TSMLT

  ! 每个应力期的输出选项
  ! IPRNT(I) = 0, 不输出到任何文件
  !          = 1, 输出到所有文件
  !          = 2, 仅输出到土壤剖面数据文件和应力期数据文件
  !          = 3, 仅输出到土壤剖面数据文件和观测点数据文件
  !          = 4, 仅输出到应力期数据文件和观测点数据文件
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPRNT

  ! 当前时段的步长
  REAL :: DELT

  ! 模拟总时长
  REAL :: TOTIM

  ! 应力期循环计数器
  INTEGER :: JPER

  ! 时段循环计数器
  INTEGER :: JSTP

  ! 迭代求解循环计数器
  INTEGER :: JITR

  ! 最大迭代次数
  INTEGER :: MAXITR

  ! 控制方程组的解是否收敛
  ! LCONV = .TURE., 收敛
  !         .FALSE., 不收敛
  LOGICAL :: LCONV

  ! 水势收敛容差 (cm)
  REAL :: HTOL

  ! 含水率收敛容差 (cm^3 cm^-3)
  REAL :: OTOL

  ! 盐分浓度收敛容差 (g cm^-3)
  REAL :: CTOL


  !------------------------------------------------------------------------------------------------
  !           定义方程组系数变量
  !------------------------------------------------------------------------------------------------
  REAL, DIMENSION(:, :), ALLOCATABLE :: CCV
  REAL, DIMENSION(:, :), ALLOCATABLE :: BBV
  REAL, DIMENSION(:, :), ALLOCATABLE :: DDV
  REAL, DIMENSION(:),    ALLOCATABLE :: COF
  REAL, DIMENSION(:),    ALLOCATABLE :: RHS


  !------------------------------------------------------------------------------------------------
  !           定义土壤材料基本特性变量
  !------------------------------------------------------------------------------------------------

  ! 干容重 (g cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: BD

  ! 砂粒占固相体积的比例
  REAL, DIMENSION(:), ALLOCATABLE :: SAND

  ! 粘粒占固相体积的比例
  REAL, DIMENSION(:), ALLOCATABLE :: CLAY

  ! 有机质占固相的比例
  REAL, DIMENSION(:), ALLOCATABLE :: OM


  !------------------------------------------------------------------------------------------------
  !           定义土壤材料水力学特性变量
  !------------------------------------------------------------------------------------------------

  ! 残余含水率 (Mualem & van Genuchten) (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: ORES

  ! 饱和含水率 (Mualem & van Genuchten) (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: OSAT

  ! 饱和导水率 (Mualem & van Genuchten) (cm d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: HYCONSAT

  ! 进气值的倒数 (Mualem & van Genuchten) (cm^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: HYPARA

  ! 孔径分布指数 (Mualem & van Genuchten)
  REAL, DIMENSION(:), ALLOCATABLE :: HYPARN

  ! 孔隙连结性参数 (Mualem & van Genuchten)
  REAL, DIMENSION(:), ALLOCATABLE :: HYPARL

  ! 凋萎含水率 (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: OWP

  ! 田间持水率 (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: OFC


  !------------------------------------------------------------------------------------------------
  !           定义土壤材料热力学特性变量
  !------------------------------------------------------------------------------------------------

  ! 砂粒的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCONSAND

  ! 粘粒的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCONCLAY

  ! 有机质的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCONOM

  ! 水的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCONWAT

  ! 空气的导热率 (J cm^-1 °C^-1 d^-1)
   REAL, DIMENSION(:), ALLOCATABLE :: THCONAIR

  ! 固相的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCONSLD

  ! 砂粒的热容量 (J cm^-3 °C^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCAPSAND

  ! 粘粒的热容量 (J cm^-3 °C^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCAPCLAY

  ! 有机质的热容量 (J cm^-3 °C^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCAPOM

  ! 水的热容量 (J cm^-3 °C^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCAPWAT

  ! 空气的热容量 (J cm^-3 °C^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCAPAIR

  ! 土壤导热率计算公式参数 (Cote & Konrad)
  REAL, DIMENSION(:), ALLOCATABLE :: THPARX

  ! 土壤导热率计算公式参数 (Cote & Konrad)
  REAL, DIMENSION(:), ALLOCATABLE :: THPARA

  ! 土壤导热率计算公式参数 (Cote & Konrad)
  REAL, DIMENSION(:), ALLOCATABLE :: THPARB


  !------------------------------------------------------------------------------------------------
  !           定义土壤材料盐分运移特性变量
  !------------------------------------------------------------------------------------------------

  ! 弥散长度 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: DISL

  ! 扩散系数 (cm^2 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: DIFW

  ! 最大吸附量 (mg kg^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: CADSMAX

  ! 吸附平衡常数 (-)
  REAL, DIMENSION(:), ALLOCATABLE :: PADS


  !------------------------------------------------------------------------------------------------
  !           定义土壤单元体基本特性变量
  !------------------------------------------------------------------------------------------------

  ! 土壤材料索引号
  INTEGER, DIMENSION(:), ALLOCATABLE :: IMAT

  ! 观测点标记
  ! IOBS(I) = 0, 不标记
  !         > 0, 用正整数标记为观测点
  INTEGER, DIMENSION(:), ALLOCATABLE :: IOBS

  ! 尺寸 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: DELV

  ! 顶部, 中心和底部的位置 (cm)
  REAL, DIMENSION(:, :), ALLOCATABLE :: LOCV


  !------------------------------------------------------------------------------------------------
  !           定义土壤单元体水分运动变量
  !------------------------------------------------------------------------------------------------

  ! 初始条件类型
  ! ISWFINI = 0, 初始条件为水势
  !         = 1, 初始条件为含水率
  INTEGER :: ISWFINI

  ! 初始水势 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: HINI

  ! 前一时段结束时的水势 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: HOLD

  ! 当前时段结束时的水势 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: HNEW

  ! 初始含水率 (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: OINI

  ! 前一时段结束时的含水率 (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: OOLD

  ! 当前时段结束时的含水率 (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: ONEW

  ! 水容量 (cm^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: HYCAP

  ! 非饱和导水率 (cm d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: HYCON

  ! 当前时段土壤单元体之间的水流量 (cm)
  REAL, DIMENSION(:, :), ALLOCATABLE :: QV

  ! 当前应力期土壤单元体的源/汇水量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WSS


  !------------------------------------------------------------------------------------------------
  !           定义土壤单元体热传导变量
  !------------------------------------------------------------------------------------------------

  ! 初始地温 (°C)
  REAL, DIMENSION(:), ALLOCATABLE :: TINI

  ! 热容量 (J cm^-3 °C^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCAP

  ! 导热率 (J cm^-1 °C^-1 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: THCON

  ! 前一时段结束时的地温 (°C)
  REAL, DIMENSION(:), ALLOCATABLE :: TOLD

  ! 当前时段结束时的地温  (°C)
  REAL, DIMENSION(:), ALLOCATABLE :: TNEW

  ! 当前时段土壤单元体之间的热流量 (J cm^-1 °C^-1)
  REAL, DIMENSION(:, :), ALLOCATABLE :: QTHV


  !------------------------------------------------------------------------------------------------
  !           定义土壤单元体盐分运移变量
  !------------------------------------------------------------------------------------------------

  ! 有效弥散系数 (cm^2 d^-1)
  REAL, DIMENSION(:, :), ALLOCATABLE :: DIS

  ! 初始盐分浓度 (mg cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: CINI
  
  ! 前一时段结束时的盐分浓度 (mg cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: COLD

  ! 当前时段结束时的盐分浓度 (mg cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: CNEW

  ! 当前时段土壤单元体之间的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:, :), ALLOCATABLE :: QSAV


  !------------------------------------------------------------------------------------------------
  !           定义土壤水分运动底部边界条件变量
  !------------------------------------------------------------------------------------------------
  ! 底部边界条件的类型
  ! ISWFBOT = 1, Dirichlet 条件, 需要输入底部边界的水势或含水率
  !         = 2, Neumann 条件, 需要输入底部边界的水量
  !         = 3, Cauchy 条件, 需要输入底部边界的水势和延迟时间
  !         = 4, 自由排水条件, 不需要输入其他数据
  !         = 5, 渗出面条件, 不需要输入其他数据
  INTEGER :: ISWFBOT

  ! 底部边界水分运动的延迟时间（d）
  INTEGER :: TLAGBOT

  ! 每个应力期底部边界的水势 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: HBOT

  ! 每个应力期底部边界的含水率 (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: OBOT

  ! 每个应力期底部边界的水流量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WBOT

  ! 当前时段底部边界的水流量 (cm)
  REAL :: QBOT


  !------------------------------------------------------------------------------------------------
  !           定义土壤热传导底部边界条件变量
  !------------------------------------------------------------------------------------------------
  
  ! 底部边界条件的类型
  ! ISHTBOT = 1, Dirichlet 条件, 需要输入底部边界的地温
  !         = 2, Neumann 条件, 需要输入底部边界的热量
  !         = 3, 零梯度条件, 不需要输入数据
  INTEGER :: ISHTBOT

  ! 底部边界的地温 (°C)
  REAL, DIMENSION(:), ALLOCATABLE :: TBOT

  ! 每个应力期底部边界的热流量 (J cm^-1 °C^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: WTHBOT

  ! 当前时段底部边界的热流量 (J cm^-1 °C^-1)
  REAL :: QTHBOT


  !------------------------------------------------------------------------------------------------
  !           定义土壤溶质运移底部边界条件变量
  !------------------------------------------------------------------------------------------------

  ! 底部边界条件的类型
  ! ISSTBOT = 1, Dirichlet 条件, 需要输入底部边界的盐分浓度
  !         = 2, Neumann 条件, 需要输入底部边界的盐分量
  !         = 3, 零梯度条件, 不需要输入数据
  INTEGER :: ISSTBOT

  ! 每个应力期底部边界的盐分浓度 (mg cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: CBOT

  ! 每个应力期底部边界的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: WSABOT

  ! 当前时段底部边界的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: QSABOT


  !------------------------------------------------------------------------------------------------
  !           定义土壤水分运动地表边界条件变量
  !------------------------------------------------------------------------------------------------
  ! 地表边界条件的类型
  ! ISWFTOP = 1, Dirichlet 条件, 需要输入地表边界的水势或含水率
  !         = 2, Neumann 条件, 需要输入地表边界的水量
  !         = 3, Cauchy 条件, 需要输入地表边界的水势和延迟使时间
  !         = 4, 大气边界
  INTEGER :: ISWFTOP

  ! 地表边界水分运动的延迟时间 (d)
  INTEGER :: TLAGTOP

  ! 每个应力期地表边界的水势 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: HTOP

  ! 每个应力期地表边界的含水率 (cm^3 cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: OTOP

  ! 每个应力期地表边界的水流量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WTOP

  ! 当前时段地表边界的水流量 (cm)
  REAL :: QTOP

  ! 地表最大允许积水深度 (cm)
  REAL :: HPONDMAX

  ! 前一时段结束时的地表积水深度 (cm)
  REAL :: HPONDOLD

  ! 当前时段结束时的地表积水深度 (cm)
  REAL :: HPONDNEW

  ! 地表风干水势 (cm)
  REAL :: HATM

  ! 每个应力期的降雨量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WRAIN

  ! 每个应力期达到冠层的水量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WCANO

  ! 每个应力期的冠层截留水量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WINT

  ! 每个应力期达到地表的水量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WSURF

  ! 每个应力期的地表径流量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WRNOF

  ! 当前时段的地表径流量 (cm)
  REAL :: QRNOF

  ! 每个应力期的地表入渗水量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WINF

  ! 当前时段的地表入渗水量 (cm)
  REAL :: QINF


  !------------------------------------------------------------------------------------------------
  !           定义土壤热传导地表边界条件变量
  !------------------------------------------------------------------------------------------------

  ! 地表边界条件的类型
  ! ISHTTOP = 1, Dirichlet 条件, 需要输入地表边界的地温
  !         = 2, Neumann 条件, 需要输入地表边界的热量
  INTEGER :: ISHTTOP

  ! 每个应力期地表边界的地温 (°C)
  REAL, DIMENSION(:), ALLOCATABLE :: TTOP

  ! 每个应力期地表边界的热流量 (J cm^-1 °C^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: WTHTOP

  ! 时段内地表边界的热流量 (J cm^-1 °C^-1)
  REAL :: QTHTOP


  !------------------------------------------------------------------------------------------------
  !           定义土壤盐分运移地表边界条件变量
  !------------------------------------------------------------------------------------------------
  ! 地表边界条件的类型
  ! ISSTTOP = 1, Dirichlet 条件, 需要输入地表边界的盐分浓度
  !         = 2, Neumann 条件, 需要输入地表边界的盐分量
  INTEGER :: ISSTTOP

  ! 每个应力期地表边界的盐分浓度 (mg cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: CTOP

  ! 每个应力期地表边界的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: WSATOP

  ! 每个应力期降雨中的盐分浓度 (mg cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: CRAIN

  ! 每个应力期降雨中的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: WSARAIN
  
  ! 每个应力期达到冠层的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: WSACANO

  ! 每个应力期冠层截留的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: WSAINT

  ! 每个应力期到达地表的水的盐分浓度 (mg cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: CSURF

  ! 每个应力期达到地表的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: WSASURF

  ! 每个应力期地表径流中的盐分流量 (mg cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: WSARNOF

  ! 当前时段地表径流中的盐分流量 (g cm^-2)
  REAL :: QSARNOF

  ! 每个应力期地表入渗的盐分流量 (g cm^-2)
  REAL, DIMENSION(:), ALLOCATABLE :: WSAINF

  ! 当前时段地表入渗的盐分流量 (g cm^-2)
  REAL :: QSAINF


  !------------------------------------------------------------------------------------------------
  !           定义大气变量
  !------------------------------------------------------------------------------------------------

  ! 太阳辐射量的计算方法
  ! IRAD = 0, 直接输入
  !      = 1, 输入日照时数计算
  INTEGER :: IRAD

  ! 纬度 (°)
  REAL :: LATITUDE

  ! 海拔 (m)
  REAL :: ALTITUDE

  ! 风速测量高度 (m)
  REAL :: WINDHGT

  ! Angstrom-Prescott 公式参数
  REAL, DIMENSION(2) :: PANGSTROM

  ! 儒略日
  INTEGER, DIMENSION(:), ALLOCATABLE :: DOY

  ! 最高气温 (°C)
  REAL, DIMENSION(:), ALLOCATABLE :: ATMAX

  ! 最低气温 (°C)
  REAL, DIMENSION(:), ALLOCATABLE :: ATMIN

  ! 平均相对湿度 (%)
  REAL, DIMENSION(:), ALLOCATABLE :: RHMEAN

  ! 平均风速 (m s^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: WSMEAN

  ! 地面 2 m 高处的平均风速 (m s^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: WS2MEAN

  ! 日照时数 (h)
  REAL, DIMENSION(:), ALLOCATABLE :: DH

  ! 最大可能日照时数 (h)
  REAL, DIMENSION(:), ALLOCATABLE :: DHMAX

  ! 太阳辐射量 (MJ m^-2 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: RAD

  ! 地外辐射量 (MJ m^-2 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: RADEXT

  ! 净辐射量 (MJ m^-2 d^-1)
  REAL, DIMENSION(:), ALLOCATABLE :: RADNET


  !------------------------------------------------------------------------------------------------
  !           定义作物冠层变量
  !------------------------------------------------------------------------------------------------

  ! 冠层消光系数
  REAL :: EXTINC

  ! 冠层反射率
  REAL :: ALBEDO

  ! 冠层截留公式经验参数 (Von Hoyningen-Hune & Braden) (cm d^-1)
  REAL :: PINT

  ! 地表覆盖度计算方法
  ! ISCF = 0, 由用户输入
  !      = 1, 输入 LAI 计算
  INTEGER :: ISCF

  ! 每个应力期的地表覆盖度
  REAL, DIMENSION(:), ALLOCATABLE :: SCF

  ! 每个应力期的叶面积指数
  REAL, DIMENSION(:), ALLOCATABLE :: LAI

  ! 每个应力期的冠层高度 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: CRPHGT

  
  !------------------------------------------------------------------------------------------------
  !           定义作物根系变量
  !------------------------------------------------------------------------------------------------

  ! 每个应力期的根系深度 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: RDEPTH

  ! 根系吸水分布的计算方法
  ! IRWUDIST = 0, 直接输入
  !          = 1, 采用指数函数计算
  INTEGER :: IRWUDIST

  ! 根系吸水分布指数函数的参数
  REAL :: PRWUDIST

  ! Feddes 模型根系吸水水分胁迫响应函数的临界水势 (cm)
  ! RWUHCRT(1) = H1,  厌氧点的临界水势, 当土壤水势大于或等于此临界值时, 水分胁迫响应函数值等于 0
  ! RWUHCRT(1) = H2,  当土壤水势介于 H3 和 H2 之间，水分胁迫响应函数值等于 1
  ! RWUHCRT(1) = H3H
  ! RWUHCRT(1) = H3L
  ! RWUHCRT(1) = H4,  凋萎时的临界水势, 当土壤水势小于或等于此临界值时, 水分胁迫响应函数值等于 0
  REAL, DIMENSION(5) :: RWUHCRT

  ! Feddes 模型根系吸水水分胁迫响应函数的临界蒸腾速率 (cm d^-1)
  ! RWUTCRT(1) = T3H
  ! RWUTCRT(2) = T3L
  REAL, DIMENSION(2) :: RWUTCRT

  ! 根系吸水补偿系数
  REAL :: RWUCOMP

  ! 根系吸水补偿系数的临界值
  REAL :: RWUCOMPCRT

  ! 当前应力期土壤剖面的根系吸水分布
  REAL, DIMENSION(:), ALLOCATABLE :: RWUDIST

  ! 当前应力期土壤剖面的根系实际吸水量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WRWUACT

  ! 当前时段土壤剖面的根系实际吸水量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: QRWUACT


  !------------------------------------------------------------------------------------------------
  !           定义蒸发蒸腾变量
  !------------------------------------------------------------------------------------------------

  ! 潜在蒸发蒸腾量的计算方法
  ! IPET = 0, 直接输入 WEPOT 和 WTPOT
  !      = 1, 输入 WETPOT, LAI 或 SCF 计算 WEPOT 和 WTPOT
  !      = 2, 采用 FAO Penman-Monteth 公式和作物系数计算 WETPOT, 再根据 LAI 或 SCF 计算 WEPOT 和 WTPOT
  !      = 3, 采用 FAO Hargreaves 公式和作物系数计算 WETPOT, 再根据 LAI 或 SCF 计算 WEPOT 和 WTPOT
  INTEGER :: IPET

  ! 每个应力期的参考蒸发蒸腾量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WETREF

  ! 每个应力期的潜在蒸发蒸腾量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WETPOT

  ! 当前时段的潜在蒸发蒸腾量 (cm)
  REAL :: QETPOT

  ! 每个应力期的地表潜在蒸发量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WEPOT

  ! 每个应力期的地表实际蒸发量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WEACT

  ! 当前时段的地表实际蒸发量 (cm)
  REAL :: QEACT

  ! 每个应力期的潜在蒸腾量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WTPOT

  ! 当前时段的潜在蒸腾量 (cm)
  REAL :: QTPOT

  ! 每个应力期的实际蒸腾量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WTACT

  ! 当前时段的实际蒸腾量 (cm)
  REAL :: QTACT

  ! 每个应力期的作物系数
  REAL, DIMENSION(:), ALLOCATABLE :: CRPCOF


  !------------------------------------------------------------------------------------------------
  !           定义灌溉变量
  !------------------------------------------------------------------------------------------------

  ! 灌溉策略
  ! IIRRSTG = 1, 固定灌溉, 根据习惯经验确定灌溉时间和水量
  !         = 2, 优化灌溉, 根据墒情指标确定灌水时间和灌水量
  INTEGER :: IIRRSTG

  ! 灌溉类别
  ! IIRRCAT = 1, 地面灌溉, 如沟灌, 喷灌, 地表滴灌等
  !         = 2, 地下灌溉, 如地下滴灌, 渗灌等
  INTEGER :: IIRRCAT

  ! 控制是否考虑作物冠层对灌溉的截留作用
  ! LIRRINT = .TRUE., 是
  !         = .FALSE., 否
  LOGICAL :: LIRRINT

  ! 灌水器所在土壤单元体的编号
  INTEGER :: ICELIRR

  ! 水分传感器所在土壤单元体的编号
  REAL :: ICELSEN

  ! 应力期内的灌水量 (cm)
  REAL, DIMENSION(:), ALLOCATABLE :: WIRR

  ! 应力期内灌溉水中的盐分浓度 (g cm^-3)
  REAL, DIMENSION(:), ALLOCATABLE :: CIRR

  ! 灌水时间的确定方法
  ! IIRRTIM = 1, 当实际蒸腾量与潜在蒸腾量的比值小于临界值时，进行灌溉
  !         = 2, 当计划湿润层土壤含水量与田间持水量的比值小于临界值时，进行灌溉
  !         = 3, 当计划湿润层土壤含水量与可利用水量的比值小于临界值时，进行灌溉
  !         = 4, 当计划湿润层土壤含水率小于临界值时，进行灌溉
  !         = 5, 当计划湿润层土壤水势小于临界值时，进行灌溉
  !         = 6, 当指定位置土壤含水量与田间持水量的比值小于临界值时，进行灌溉
  !         = 7, 当指定位置土壤含水量与可利用水量的比值小于临界值时，进行灌溉
  !         = 8, 当指定位置土壤含水率小于临界值时，进行灌溉
  !         = 9, 当指定位置土壤水势小于临界值时，进行灌溉
  INTEGER :: IIRRTIM

  ! 确定灌水时间的指标临界值
  ! 当 IIRRTIM = 1 时, IRRTIMCRT 为实际蒸腾量与潜在蒸腾量比值的临界值
  !            = 2 时, IRRTIMCRT 为计划湿润层土壤含水量与田间持水量比值的临界值
  !            = 3 时, IRRTIMCRT 为计划湿润层土壤含水量与可利用水量比值的临界值
  !            = 4 时, IRRTIMCRT 为计划湿润层土壤含水率的临界值
  !            = 5 时, IRRTIMCRT 为计划湿润层土壤水势的临界值
  !            = 6 时, IRRTIMCRT 为指定位置传感器土壤含水量与田间持水量比值的临界值
  !            = 7 时, IRRTIMCRT 为指定位置传感器土壤含水量与可利用水量比值的临界值
  !            = 8 时, IRRTIMCRT 为指定位置传感器土壤含水率的临界值
  !            = 9 时, IRRTIMCRT 为指定位置传感器土壤水势的临界值
  REAL :: IRRTIMCRT

  ! 灌水量的确定方法
  ! IIRRVOL = 1, 由用户直接输入灌水量
  !         = 2, 根据计划湿润层田间持水量的比例确定
  !         = 3, 根据潜在蒸发蒸腾量的比例确定
  INTEGER :: IIRRVOL

  ! 确定灌水量的指标临界值
  ! 当 IIRRVOL = 1 时, IRRVOLCRT 为灌水量, cm
  !            = 2 时, IRRVOLCRT 为计划湿润层目标土壤含水量与田间持水量比值的临界值
  !            = 3 时, IRRVOLCRT 为目标灌水量与潜在蒸发蒸腾量比值的临界值
  REAL :: IRRVOLCRT

  ! 计划湿润层深度 (cm)
  REAL :: IRRDEPTH

  ! 控制是否考虑土壤固相对溶质的吸附作用
  ! LADS = .TRUE., 是
  !      = .FALSE., 否
  LOGICAL :: LADS

  ! 弥散系数弯曲因子的计算方法
  ! ITOR = 1, 采用 Modrup 公式计算
  !      = 2, 采用 Millington & Quirk 公式计算
  INTEGER :: ITOR

  ! 时段内的降雨量 (cm)
  REAL :: QRAIN

  ! 时段内达到冠层的水量 (cm)
  REAL :: QCANO

  ! 时段内的作物冠层截留量 (cm)
  REAL :: QINT


  ! 时段内达到地表的水量 (cm)
  REAL :: QSURF


  ! 时段内的地表潜在蒸发量 (cm)
  REAL :: QEPOT

END MODULE VAR

