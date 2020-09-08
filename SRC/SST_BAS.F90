!==================================================================================================
! 子程序包: SST_BAS (Basic Package for Soil Salt Transport Process)
! 主要功能: 计算土壤盐分运移参数
! 创建日期: 2020年7月4日
! 修改日志:
!==================================================================================================


SUBROUTINE SST_BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('ISSTBOT')) THEN
    CALL RDSINT('ISSTBOT', ISSTBOT)
  ELSE
    ISSTBOT = 0
  END IF
  IF (RDINQR('ISSTTOP')) THEN
    CALL RDSINT('ISSTTOP', ISSTTOP)
  ELSE
    ISSTTOP = 0
  END IF
  IF (RDINQR('ITOR')) THEN
    CALL RDSINT('ITOR', ITOR)
  ELSE
    ITOR = 1
  END IF

  L_SST_STO = .TRUE.
  L_SST_CCD = .TRUE.
  L_SST_CCA = .TRUE.
  L_SST_GBB = .TRUE.
  L_SST_GTB = (ISWFTOP /= 4)
  L_SST_INF = (ISWFTOP == 4)
  L_SST_TMA = .TRUE.

  RETURN
END SUBROUTINE SST_BAS_DF


SUBROUTINE SST_BAS_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(BD))      ALLOCATE(BD(NMAT),      SOURCE = 0.0)
  IF (.NOT. ALLOCATED(DISL))    ALLOCATE(DISL(NMAT),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(DIFW))    ALLOCATE(DIFW(NMAT),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(DIS))     ALLOCATE(DIS(NCEL, 2),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(COF))     ALLOCATE(COF(NCEL),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RHS))     ALLOCATE(RHS(NCEL),     SOURCE = 0.0)

  RETURN
END SUBROUTINE SST_BAS_AL


SUBROUTINE SST_BAS_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 局部变量
  INTEGER :: I

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('BD')) THEN
    CALL RDFREA('BD', BD, NMAT, NMAT)
  ELSE
    BD(:) = 1.5
  END IF
  IF (RDINQR('DISL')) THEN
    CALL RDFREA('DISL', DISL, NMAT, NMAT)
  ELSE
    DISL(:) = 10.0
  END IF
  IF (RDINQR('DIFW')) THEN
    CALL RDFREA('DIFW', DIFW, NMAT, NMAT)
  ELSE
    DIFW(:) = 1.0
  END IF

  RETURN
END SUBROUTINE SST_BAS_RP


SUBROUTINE SST_BAS_FM
!**************************************************************************************************
! 计算盐分运移参数和方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J, J1, J2
  REAL :: O1, O2
  REAL :: OSAT1, OSAT2
  REAL :: DISL1, DISL2
  REAL :: DIFW1, DIFW2
  REAL :: DELV1, DELV2
  REAL :: TOR1, TOR2        ! Tortuosity factor

  ! 第 2 至 NCEL - 1 单元体
  DO I = 2, NCEL - 1
    J = IMAT(I)
    J1 = IMAT(I - 1)
    J2 = IMAT(I + 1)
    O1 = (DELV(I - 1) * ONEW(I - 1) + DELV(I) * ONEW(I)) / (DELV(I - 1) + DELV(I))
    O2 = (DELV(I + 1) * ONEW(I + 1) + DELV(I) * ONEW(I)) / (DELV(I + 1) + DELV(I))
    OSAT1 = (DELV(I - 1) * OSAT(J1) + DELV(I) * OSAT(J)) / (DELV(I - 1) + DELV(I))
    OSAT2 = (DELV(I + 1) * OSAT(J2) + DELV(I) * OSAT(J)) / (DELV(I + 1) + DELV(I))
    DISL1 = (DELV(I - 1) * DISL(J1) + DELV(I) * DISL(J)) / (DELV(I - 1) + DELV(I))
    DISL2 = (DELV(I + 1) * DISL(J2) + DELV(I) * DISL(J)) / (DELV(I + 1) + DELV(I))
    DIFW1 = (DELV(I - 1) * DIFW(J1) + DELV(I) * DIFW(J)) / (DELV(I - 1) + DELV(I))
    DIFW2 = (DELV(I + 1) * DIFW(J2) + DELV(I) * DIFW(J)) / (DELV(I + 1) + DELV(I))

    SELECT CASE (ITOR)
    CASE DEFAULT
      TOR1 = 0.66 * (O1 / OSAT1) ** (8.0 / 3.0)
      TOR2 = 0.66 * (O2 / OSAT2) ** (8.0 / 3.0)
    CASE (2)
      TOR1 = O1 ** (7.0 / 3.0) / OSAT1 ** 2.0
      TOR2 = O2 ** (7.0 / 3.0) / OSAT2 ** 2.0
    END SELECT
    
    DIS(I, 1) = DISL1 * ABS(QV(I, 1) / DELT / O1) + O1 * DIFW1 * TOR1
    DIS(I, 2) = DISL2 * ABS(QV(I, 2) / DELT / O2) + O2 * DIFW2 * TOR2

    COF(I) = 0.0
    RHS(I) = 0.0
  END DO

  ! 第 1 单元体
  I = 1
  J = IMAT(I)
  J2 = IMAT(I + 1)
  O2 = (DELV(I + 1) * ONEW(I + 1) + DELV(I) * ONEW(I)) / (DELV(I + 1) + DELV(I))
  OSAT2 = (DELV(I + 1) * OSAT(J2) + DELV(I) * OSAT(J)) / (DELV(I + 1) + DELV(I))
  DISL2 = (DELV(I + 1) * DISL(J2) + DELV(I) * DISL(J)) / (DELV(I + 1) + DELV(I))
  DIFW2 = (DELV(I + 1) * DIFW(J2) + DELV(I) * DIFW(J)) / (DELV(I + 1) + DELV(I))
  
  SELECT CASE (ITOR)
  CASE DEFAULT
    TOR2 = 0.66 * (O2 / OSAT2) ** (8.0 / 3.0)
  CASE (2)
    TOR2 = O2 ** (7.0 / 3.0) / OSAT2 ** 2.0
  END SELECT
  
  DIS(I, 2) = DISL2 * ABS(QV(I, 2) / DELT / O2) + O2 * DIFW2 * TOR2
  
  COF(I) = 0.0
  RHS(I) = 0.0

  ! 第 NCEL 单元体
  I = NCEL
  J1 = IMAT(I - 1)
  O1 = (DELV(I - 1) * ONEW(I - 1) + DELV(I) * ONEW(I)) / (DELV(I - 1) + DELV(I))
  OSAT1 = (DELV(I - 1) * OSAT(J1) + DELV(I) * OSAT(J)) / (DELV(I - 1) + DELV(I))
  DISL1 = (DELV(I - 1) * DISL(J1) + DELV(I) * DISL(J)) / (DELV(I - 1) + DELV(I))
  DIFW1 = (DELV(I - 1) * DIFW(J1) + DELV(I) * DIFW(J)) / (DELV(I - 1) + DELV(I))

  SELECT CASE (ITOR)
  CASE DEFAULT
    TOR1 = 0.66 * (O1 / OSAT1) ** (8.0 / 3.0)
  CASE (2)
    TOR1 = O1 ** (7.0 / 3.0) / OSAT1 ** 2.0
  END SELECT

  DIS(I, 1) = DISL1 * ABS(QV(I, 1) / DELT / O1) + O1 * DIFW1 * TOR1

  COF(I) = 0.0
  RHS(I) = 0.0

  RETURN
END SUBROUTINE SST_BAS_FM


SUBROUTINE SST_BAS_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(BD))      DEALLOCATE(BD)
  IF (ALLOCATED(DISL))    DEALLOCATE(DISL)
  IF (ALLOCATED(DIFW))    DEALLOCATE(DIFW)
  IF (ALLOCATED(DIS))     DEALLOCATE(DIS)
  IF (ALLOCATED(COF))     DEALLOCATE(COF)
  IF (ALLOCATED(RHS))     DEALLOCATE(RHS)

  RETURN
END SUBROUTINE SST_BAS_DA

