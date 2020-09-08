!==================================================================================================
! 子程序包: IRR_OWA (Optimal Water Application Package for Irrigation Process)
! 主要功能: 模拟优化灌溉制度
! 创建日期: 2020年7月14日
! 修改日志:
!==================================================================================================


SUBROUTINE IRR_OWA_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('ICELSEN')) THEN
    CALL RDSINT('ICELSEN', ICELSEN)
  ELSE
    ICELSEN = 1
  END IF
  IF (RDINQR('IIRRTIM')) THEN
    CALL RDSINT('IIRRTIM', IIRRTIM)
  ELSE
    IIRRTIM = 3
  END IF
  IF (RDINQR('IRRTIMCRT')) THEN
    CALL RDSREA('IRRTIMCRT', IRRTIMCRT)
  ELSE
    IRRTIMCRT = 1.0
  END IF
  IF (RDINQR('IIRRVOL')) THEN
    CALL RDSINT('IIRRVOL', IIRRVOL)
  ELSE
    IIRRVOL = 2
  END IF
  IF (RDINQR('IRRVOLCRT')) THEN
    CALL RDSINT('IRRVOLCRT', IRRVOLCRT)
  ELSE
    IRRVOLCRT = 1.0
  END IF

  RETURN
END SUBROUTINE IRR_OWA_DF


SUBROUTINE IRR_OWA_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE IRR_OWA_AL


SUBROUTINE IRR_OWA_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR
  
  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('WIRR')) THEN
    CALL RDFREA('WIRR', WIRR, NPER, NPER)
  ELSE
    WIRR(:) = 0.0
  END IF  

  RETURN
END SUBROUTINE IRR_OWA_RP


SUBROUTINE IRR_OWA_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  LOGICAL :: LDOIRR
  INTEGER :: I, J
  REAL    :: SUM1, SUM2

  LDOIRR = .FALSE.

  SELECT CASE (IIRRTIM)
  
  ! 当实际蒸腾量与潜在蒸腾量的比值小于临界值时，进行灌溉
  CASE (1)
    IF (WTACT(JPER) < WTPOT(JPER) * IRRTIMCRT) LDOIRR = .TRUE.
  
  ! 当计划湿润层土壤含水量与田间持水量的比值小于临界值时，进行灌溉
  CASE (2) 
    SUM1 = 0.0
    SUM2 = 0.0
    DO I = 1, NCEL
      IF (LOCV(I, 2) <= IRRDEPTH) THEN
        J = IMAT(I)
        SUM1 = SUM1 + ONEW(I) * DELV(I)
        SUM2 = SUM2 + OFC(J) * DELV(I) * IRRTIMCRT
      END IF
    END DO
    IF (SUM1 <= SUM2) LDOIRR = .TRUE.

  ! 当计划湿润层土壤含水量与可利用水量的比值小于临界值时，进行灌溉
  CASE (3)
    SUM1 = 0.0
    SUM2 = 0.0
    DO I = 1, NCEL
      IF (LOCV(I, 2) <= IRRDEPTH) THEN
        J = IMAT(I)
        SUM1 = SUM1 + ONEW(I) * DELV(I)
        SUM2 = SUM2 + (OFC(J) - OWP(J)) * DELV(I) * IRRTIMCRT
      END IF
    END DO
    IF (SUM1 <= SUM2) LDOIRR = .TRUE.

  ! 当计划湿润层土壤含水率小于临界值时，进行灌溉
  CASE (4)
    SUM1 = 0.0
    SUM2 = 0.0
    DO I = 1, NCEL
      IF (LOCV(I, 2) <= IRRDEPTH) THEN
        SUM1 = SUM1 + ONEW(I) * DELV(I)
        SUM2 = SUM2 + IRRTIMCRT * DELV(I)
      END IF
    END DO
    IF (SUM1 < SUM2) LDOIRR = .TRUE.

  ! 当计划湿润层土壤水势小于临界值时，进行灌溉
  CASE (5)
    SUM1 = 0.0
    SUM2 = 0.0
    DO I = 1, NCEL
      IF (LOCV(I, 2) <= IRRDEPTH) THEN
        J = IMAT(I)
        SUM1 = SUM1 + HNEW(I) * DELV(I)
        SUM2 = SUM2 + IRRTIMCRT * DELV(I)
      END IF
    END DO
    IF (SUM1 <= SUM2) LDOIRR = .TRUE.

  ! 当指定位置土壤含水量与田间持水量的比值小于临界值时，进行灌溉
  CASE (6)
    J = IMAT(ICELSEN)
    IF (ONEW(ICELSEN) <= OFC(J) * IRRTIMCRT) LDOIRR = .TRUE.

  ! 当指定位置土壤含水量与可利用水量的比值小于临界值时，进行灌溉
  CASE (7)
    J = IMAT(ICELSEN)
    IF (HNEW(ICELSEN) <= (OFC(J) - OWP(J)) * IRRTIMCRT) LDOIRR = .TRUE.

  ! 当指定位置土壤含水率小于临界值时，进行灌溉
  CASE (8)
    IF (ONEW(ICELSEN) <= IRRTIMCRT) LDOIRR = .TRUE.

  ! 当指定位置土壤水势小于临界值时，进行灌溉
  CASE (9)
    IF (HNEW(ICELSEN) <= IRRTIMCRT) LDOIRR = .TRUE.

  END SELECT

  IF (LDOIRR) THEN

    SELECT CASE (IIRRVOL)

    ! 由用户直接输入灌水量
    CASE (1)

    ! 根据计划湿润层田间持水量的比例确定
    CASE (2)
      SUM1 = 0.0
      DO I = 1, NCEL
        IF (IRRDEPTH >= LOCV(I, 2)) THEN
          J = IMAT(I)
          SUM1 = SUM1 + (OFC(J) * IRRTIMCRT - ONEW(I)) * DELV(I)
        END IF
      END DO
      WIRR(JPER) = SUM1

    ! 根据潜在蒸发蒸腾量的比例确定
    CASE (3)
      WIRR(JPER) = WETPOT(JPER) * IRRVOLCRT
    END SELECT

    IF (IIRRSTG == 1) THEN
      IF (LIRRINT) THEN
        ! 考虑作物冠层对灌溉的截留作用
        WCANO(JPER) = WCANO(JPER) + WIRR(JPER)
      ELSE
        ! 不考虑作物冠层对灌溉的截留作用
        WSURF(JPER) = WSURF(JPER) + WIRR(JPER)
      END IF
    END IF

  END IF

  RETURN
END SUBROUTINE IRR_OWA_ST


SUBROUTINE IRR_OWA_FM
!**************************************************************************************************
! 计算方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (IIRRCAT == 2) THEN
    IF ((ICELIRR > 0) .AND. (ICELIRR <= NCEL)) &
      RHS(ICELIRR) = RHS(ICELIRR) + DELT * WIRR(JPER)
  END IF

  RETURN
END SUBROUTINE IRR_OWA_FM


SUBROUTINE IRR_OWA_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE IRR_OWA_DA

