!==================================================================================================
! 子程序包: SWF_BAS (Basic Package for Soil Water Flow Process)
! 主要功能: 设置土壤水分运动的边界条件类型, 计算土壤水分运动参数
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE SWF_BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('ISWFBOT')) THEN
    CALL RDSINT('ISWFBOT', ISWFBOT)
  ELSE
    ISWFBOT = 0
  END IF
  IF (RDINQR('ISWFTOP')) THEN
    CALL RDSINT('ISWFTOP', ISWFTOP)
  ELSE
    ISWFTOP = 0
  END IF

  ! 设置子程序包
  L_SWF_STO = .TRUE.
  L_SWF_CCF = .TRUE.
  L_SWF_GBB = .TRUE.
  L_SWF_GTB = (ISWFTOP /= 4)
  L_SWF_INF = (ISWFTOP == 4)
  L_SWF_EVP = (ISWFTOP == 4)
  L_SWF_RWU = LRWU
  L_SWF_TMA = .TRUE.

  RETURN
END SUBROUTINE SWF_BAS_DF


SUBROUTINE SWF_BAS_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(ORES))     ALLOCATE(ORES(NMAT),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OSAT))     ALLOCATE(OSAT(NMAT),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HYCONSAT)) ALLOCATE(HYCONSAT(NMAT), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HYPARA))   ALLOCATE(HYPARA(NMAT),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HYPARN))   ALLOCATE(HYPARN(NMAT),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HYPARL))   ALLOCATE(HYPARL(NMAT),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HYCAP))    ALLOCATE(HYCAP(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HYCON))    ALLOCATE(HYCON(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(COF))      ALLOCATE(COF(NCEL),      SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RHS))      ALLOCATE(RHS(NCEL),      SOURCE = 0.0)


  RETURN
END SUBROUTINE SWF_BAS_AL


SUBROUTINE SWF_BAS_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 局部变量
  INTEGER :: I, J

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('ORES')) THEN
    CALL RDFREA('ORES', ORES, NMAT, NMAT)
  ELSE
    ORES(:) = 0.078
  END IF
  IF (RDINQR('OSAT')) THEN
    CALL RDFREA('OSAT', OSAT, NMAT, NMAT)
  ELSE
    OSAT(:) = 0.43
  END IF
  IF (RDINQR('HYCONSAT')) THEN
    CALL RDFREA('HYCONSAT', HYCONSAT, NMAT, NMAT)
  ELSE
    HYCONSAT(:) = 24.9
  END IF
  IF (RDINQR('HYPARA')) THEN
    CALL RDFREA('HYPARA', HYPARA, NMAT, NMAT)
  ELSE
    HYPARA(:) = 0.036
  END IF
  IF (RDINQR('HYPARN')) THEN
    CALL RDFREA('HYPARN', HYPARN, NMAT, NMAT)
  ELSE
    HYPARN(:) = 1.56
  END IF
  IF (RDINQR('HYPARL')) THEN
    CALL RDFREA('HYPARL', HYPARL, NMAT, NMAT)
  ELSE
    HYPARL(:) = 0.5
  END IF

  RETURN
END SUBROUTINE SWF_BAS_RP


SUBROUTINE SWF_BAS_FM
!**************************************************************************************************
! 计算土壤水分运动参数并重置 COF 和 RHS
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J

  DO I = 1, NCEL
    J = IMAT(I)
    CALL U_O_MVG(ONEW(I), HNEW(I), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    CALL U_HYCAP_MVG(HYCAP(I), HNEW(I), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    CALL U_HYCON_MVG(HYCON(I), HNEW(I), HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
    COF(I) = 0.0
    RHS(I) = 0.0
  END DO

  RETURN
END SUBROUTINE SWF_BAS_FM


SUBROUTINE SWF_BAS_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(ORES))     DEALLOCATE(ORES)
  IF (ALLOCATED(OSAT))     DEALLOCATE(OSAT)
  IF (ALLOCATED(HYCONSAT)) DEALLOCATE(HYCONSAT)
  IF (ALLOCATED(HYPARA))   DEALLOCATE(HYPARA)
  IF (ALLOCATED(HYPARN))   DEALLOCATE(HYPARN)
  IF (ALLOCATED(HYPARL))   DEALLOCATE(HYPARL)
  IF (ALLOCATED(HYCAP))    DEALLOCATE(HYCAP)
  IF (ALLOCATED(HYCON))    DEALLOCATE(HYCON)
  IF (ALLOCATED(COF))      DEALLOCATE(COF)
  IF (ALLOCATED(RHS))      DEALLOCATE(RHS)

  RETURN
END SUBROUTINE SWF_BAS_DA

