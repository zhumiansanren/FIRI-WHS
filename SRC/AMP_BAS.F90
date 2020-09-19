!==================================================================================================
! 子程序包: AMP_BAS (Basic Package for Atmosphere Process)
! 主要功能: 为模型提供纬度, 海拔, 气温, 风速, 湿度, 辐射等数据
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE AMP_BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('LATITUDE')) THEN
    CALL RDSREA('LATITUDE', LATITUDE)
  ELSE
    LATITUDE = 35.0
  END IF
  IF (RDINQR('ALTITUDE')) THEN
    CALL RDSREA('ALTITUDE', ALTITUDE)
  ELSE
    ALTITUDE = 100.0
  END IF
  IF (RDINQR('WINDHGT')) THEN
    CALL RDSREA('WINDHGT', WINDHGT)
  ELSE
    WINDHGT = 10.0
  END IF
  IF (RDINQR('PANGSTROM')) THEN
    CALL RDFREA('PANGSTROM', PANGSTROM, 2, 2)
  ELSE
    PANGSTROM(:) = (/0.25, 0.5/)
  END IF
  IF (RDINQR('IRAD')) THEN
    CALL RDSINT('IRAD', IRAD)
  ELSE
    IRAD = 1
  END IF

  L_AMP_PCP = .TRUE.
  L_AMP_PET = .TRUE.

  RETURN
END SUBROUTINE AMP_BAS_DF


SUBROUTINE AMP_BAS_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(DOY))     ALLOCATE(DOY(NPER),     SOURCE = 0)
  IF (.NOT. ALLOCATED(ATMAX))   ALLOCATE(ATMAX(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(ATMIN))   ALLOCATE(ATMIN(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RHMEAN))  ALLOCATE(RHMEAN(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSMEAN))  ALLOCATE(WSMEAN(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WS2MEAN)) ALLOCATE(WS2MEAN(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(DH))      ALLOCATE(DH(NPER),      SOURCE = 0.0)
  IF (.NOT. ALLOCATED(DHMAX))   ALLOCATE(DHMAX(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RAD))     ALLOCATE(RAD(NPER),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RADEXT))  ALLOCATE(RADEXT(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RADNET))  ALLOCATE(RADNET(NPER),  SOURCE = 0.0)

  RETURN
END SUBROUTINE AMP_BAS_AL


SUBROUTINE AMP_BAS_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_MET 文件读取数据
  CALL RDINIT(U_IN_MET, U_LOG, S_IN_MET)
  IF (RDINQR('DOY')) THEN
    CALL RDFINT('DOY', DOY, NPER, NPER)
  ELSE
    DOY(:) = 100
  END IF
  IF (RDINQR('ATMAX')) THEN
    CALL RDFREA('ATMAX', ATMAX, NPER, NPER)
  ELSE
    ATMAX(:) = 25.0
  END IF
  IF (RDINQR('ATMIN')) THEN
    CALL RDFREA('ATMIN', ATMIN, NPER, NPER)
  ELSE
    ATMIN(:) = 15.0
  END IF
  IF (RDINQR('RHMEAN')) THEN
    CALL RDFREA('RHMEAN', RHMEAN, NPER, NPER)
  ELSE
    RHMEAN(:) = 75.0
  END IF
  IF (RDINQR('WSMEAN')) THEN
    CALL RDFREA('WSMEAN', WSMEAN, NPER, NPER)
  ELSE
    WSMEAN(:) = 2.0
  END IF
  IF (RDINQR('DH')) THEN
    CALL RDFREA('DH', DH, NPER, NPER)
  ELSE
    DH(:) = 8.0
  END IF
  IF (RDINQR('RAD')) THEN
    CALL RDFREA('RAD', RAD, NPER, NPER)
  ELSE
    RAD(:) = 15.0
  END IF

  RETURN
END SUBROUTINE AMP_BAS_RP


SUBROUTINE AMP_BAS_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 常数:
  REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)   ! 圆周率
  REAL, PARAMETER :: SB = 4.903E-9          ! Stefan-Boltzmann 常数 (MJ K^-4 m^-2 d^-1)

  ! 函数声明:
  REAL :: SVP

  ! 局部变量
  REAL :: Y                                 ! 纬度, 以弧度表示 (rad)
  REAL :: DR                                ! 相对日地距离的倒数
  REAL :: D                                 ! 太阳偏角 (rad)
  REAL :: OMG                               ! 日落时角 (rad)
  REAL :: RSO                               ! 晴空辐射量 (MJ m^-2 d^-1)
  REAL :: RNS                               ! 净短波辐射量 (MJ m^-2 d^-1)
  REAL :: RNL                               ! 净长波辐射量 (MJ m^-2 d^-1)
  REAL :: ES                                ! 饱和水汽压 (kPa)
  REAL :: EA                                ! 实际水汽压 (kPa)

  ! 计算地外辐射量和最大可能日照时数
  Y = LATITUDE * PI / 180.0 ! 角度转化为弧度
  DR = 1.0 + 0.033 * COS(2.0 * PI * DOY(JPER) / 365.0)
  D = 0.409 * SIN(2.0 * PI * DOY(JPER) / 365.0 - 1.39)
  OMG = ACOS(-TAN(Y) * TAN(D))
  RADEXT(JPER) = 118.08 / PI * DR * (OMG * SIN(Y) * SIN(D) + SIN(OMG) * COS(Y) * COS(D))
  DHMAX(JPER) = 24.0 / PI * OMG

  ! 如果输入的是日照时数, 则需要计算太阳辐射量
  IF (IRAD == 1) THEN
    RAD(JPER) = (PANGSTROM(1) + PANGSTROM(2) * DH(JPER) / DHMAX(JPER)) * RADEXT(JPER)
  END IF

  ! 计算净辐射量
  ES = (SVP(ATMAX(JPER)) + SVP(ATMIN(JPER))) * 0.5
  EA = ES * RHMEAN(JPER) * 0.01
  RSO = (PANGSTROM(1) + PANGSTROM(2)) * RADEXT(JPER)
  RNS = (1.0 - ALBEDO) * RAD(JPER)
  RNL = 4.903E-9 * ((ATMAX(JPER) + 273.16) ** 4.0 + (ATMIN(JPER) + 273.16) ** 4.0) * 0.5 * &
       (0.34 - 0.14 * SQRT(EA)) * (1.35 * RAD(JPER) / RSO - 0.35)
  RADNET(JPER) = RNS - RNL

  ! 计算地面 2 m 高处的风速
  IF (WINDHGT > 2.0) THEN
    WS2MEAN(JPER) = WSMEAN(JPER) * 4.87 / LOG(67.8 * WINDHGT - 5.42)
  END IF

  RETURN
END SUBROUTINE AMP_BAS_ST


SUBROUTINE AMP_BAS_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(DOY))     DEALLOCATE(DOY)
  IF (ALLOCATED(ATMAX))   DEALLOCATE(ATMAX)
  IF (ALLOCATED(ATMIN))   DEALLOCATE(ATMIN)
  IF (ALLOCATED(RHMEAN))  DEALLOCATE(RHMEAN)
  IF (ALLOCATED(WSMEAN))  DEALLOCATE(WSMEAN)
  IF (ALLOCATED(WS2MEAN)) DEALLOCATE(WS2MEAN)
  IF (ALLOCATED(DH))      DEALLOCATE(DH)
  IF (ALLOCATED(DHMAX))   DEALLOCATE(DHMAX)
  IF (ALLOCATED(RAD))     DEALLOCATE(RAD)
  IF (ALLOCATED(RADEXT))  DEALLOCATE(RADEXT)
  IF (ALLOCATED(RADNET))  DEALLOCATE(RADNET)

  RETURN
END SUBROUTINE AMP_BAS_DA

