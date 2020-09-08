!==================================================================================================
! 子程序包: CRP_PET (Potential Evaportranspiration Package for Crop Process)
! 主要功能: 计算潜在蒸发蒸腾量
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE CRP_PET_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('IPET')) THEN
    CALL RDSINT('IPET', IPET)
  ELSE
    IPET = 0
  END IF

  RETURN
END SUBROUTINE CRP_PET_DF


SUBROUTINE CRP_PET_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(SCF))    ALLOCATE(SCF(NPER),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(CRPCOF)) ALLOCATE(CRPCOF(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WETPOT)) ALLOCATE(WETPOT(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WEPOT))  ALLOCATE(WEPOT(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WTPOT))  ALLOCATE(WTPOT(NPER),  SOURCE = 0.0)

  RETURN
END SUBROUTINE CRP_PET_AL


SUBROUTINE CRP_PET_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_MET 文件读取数据
  CALL RDINIT(U_IN_MET, U_LOG, S_IN_MET)
  IF (RDINQR('CRPCOF')) THEN
    CALL RDFREA('CRPCOF', CRPCOF, NPER, NPER)
  ELSE
    CRPCOF(:) = 1.0
  END IF
  IF (RDINQR('WETPOT')) THEN
    CALL RDFREA('WETPOT', WETPOT, NPER, NPER)
  END IF

  RETURN
END SUBROUTINE CRP_PET_RP


SUBROUTINE CRP_PET_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  REAL :: SVP

  ! 局部变量
  REAL :: ATMEAN             ! 平均气温 (°C)
  REAL :: GAMMA              ! Psychrometric constant (kPa °C^-1)
  REAL :: DELTA              ! Slope vapour pressure curve (kPa °C^-1)
  REAL :: ES                 ! 饱和水汽压 (kPa)
  REAL :: EA                 ! 实际水汽压 (kPa)
  REAL :: ETRAD              ! 由太阳辐射产生的潜在蒸腾蒸发速率 (mm d^-1)
  REAL :: ETAERO             ! 由空气流动产生的潜在蒸腾蒸发速率 (mm d^-1)
  REAL :: ET0                ! 参考作物蒸腾蒸发速率 (mm d^-1)

  
  SELECT CASE (IPET)

  CASE(1)

    WEPOT(JPER) = WETPOT(JPER) * (1.0 - SCF(JPER))
    WTPOT(JPER) = WETPOT(JPER) - WEPOT(JPER)

  ! 采用 FAO Penman-Monteith 公式计算参考作物蒸发蒸腾速率
  CASE (2)
    ATMEAN = (ATMIN(JPER) + ATMAX(JPER)) * 0.5
    GAMMA = 0.000665 * 101.3 * (1.0 - 0.0065 * ALTITUDE / 293.0) ** 5.26
    DELTA = 4098.0 * SVP(ATMEAN) / (ATMEAN + 237.3) ** 2.0
    ES = (SVP(ATMIN(JPER)) + SVP(ATMAX(JPER))) * 0.5
    EA = ES * RHMEAN(JPER) * 0.01
    ETRAD = 0.408 * DELTA * RADNET(JPER) / (DELTA + GAMMA * (1.0 + 0.34 * WS2MEAN(JPER)))
    ETAERO = GAMMA * 900.0 / (ATMEAN + 273.0) * WS2MEAN(JPER) * (ES - EA) / &
            (DELTA + GAMMA * (1.0 + 0.34 * WS2MEAN(JPER)))
    ET0 = ETRAD + ETAERO
    IF (ET0 < 0.0) ET0 = 0.0

    WETPOT(JPER) = ET0 * 0.1 * CRPCOF(JPER)
    WEPOT(JPER) = WETPOT(JPER) * (1.0 - SCF(JPER))
    WTPOT(JPER) = WETPOT(JPER) - WEPOT(JPER)

  ! 采用 FAO Hargreaves 公式计算参考蒸发蒸腾速率
  CASE (3)
    ATMEAN = (ATMIN(JPER) + ATMAX(JPER)) * 0.5
    ET0 = 0.0023 * (ATMEAN + 17.8) * SQRT(ATMAX(JPER) - ATMIN(JPER)) * RADEXT(JPER)
    IF (ET0 < 0.0) ET0 = 0.0

    WETPOT(JPER) = ET0 * 0.1 * CRPCOF(JPER)
    WEPOT(JPER) = WETPOT(JPER) * (1.0 - SCF(JPER))
    WTPOT(JPER) = WETPOT(JPER) - WEPOT(JPER)

  END SELECT

  RETURN
END SUBROUTINE CRP_PET_ST


SUBROUTINE CRP_PET_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(SCF))    DEALLOCATE(SCF)
  IF (ALLOCATED(CRPCOF)) DEALLOCATE(CRPCOF)
  IF (ALLOCATED(WETPOT)) DEALLOCATE(WETPOT)
  IF (ALLOCATED(WEPOT))  DEALLOCATE(WEPOT)
  IF (ALLOCATED(WTPOT))  DEALLOCATE(WTPOT)

  RETURN
END SUBROUTINE CRP_PET_DA

