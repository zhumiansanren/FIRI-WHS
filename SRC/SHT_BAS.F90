!==================================================================================================
! 子程序包: SHT_BAS (Basic Package for Soil Heat Transport Process)
! 主要功能: 设置土壤热传导的边界条件类型, 计算土壤热传导参数
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE SHT_BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('ISHTBOT')) THEN
    CALL RDSINT('ISHTBOT', ISHTBOT)
  ELSE
    ISHTBOT = 0
  END IF
  IF (RDINQR('ISHTTOP')) THEN
    CALL RDSINT('ISHTTOP', ISHTTOP)
  ELSE
    ISHTTOP = 0
  END IF

  L_SHT_STO = .TRUE.
  L_SHT_CCT = .TRUE.
  L_SHT_GBB = .TRUE.
  L_SHT_GTB = .TRUE.
  L_SHT_TMA = .TRUE.

  RETURN
END SUBROUTINE SHT_BAS_DF


SUBROUTINE SHT_BAS_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(BD))       ALLOCATE(BD(NMAT),       SOURCE = 0.0)
  IF (.NOT. ALLOCATED(SAND))     ALLOCATE(SAND(NMAT),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(CLAY))     ALLOCATE(CLAY(NMAT),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(THCONSLD)) ALLOCATE(THCONSLD(NMAT), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(THCONWAT)) ALLOCATE(THCONWAT(NMAT), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(THPARX))   ALLOCATE(THPARX(NMAT),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(THPARA))   ALLOCATE(THPARA(NMAT),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(THPARB))   ALLOCATE(THPARB(NMAT),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(THCON))    ALLOCATE(THCON(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(THCAP))    ALLOCATE(THCAP(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(COF))      ALLOCATE(COF(NCEL),      SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RHS))      ALLOCATE(RHS(NCEL),      SOURCE = 0.0)

  RETURN
END SUBROUTINE SHT_BAS_AL


SUBROUTINE SHT_BAS_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('BD')) THEN
    CALL RDFREA('BD', BD, NMAT, NMAT)
  ELSE
    BD(:) = 1.5
  END IF
  IF (RDINQR('SAND')) THEN
    CALL RDFREA('SAND', SAND, NMAT, NMAT)
  ELSE
    SAND(:) = 0.5
  END IF
  IF (RDINQR('CLAY')) THEN
    CALL RDFREA('CLAY', CLAY, NMAT, NMAT)
  ELSE
    CLAY(:) = 0.2
  END IF
  IF (RDINQR('THCONSLD')) THEN
    CALL RDFREA('THCONSLD', THCONSLD, NMAT, NMAT)
  ELSE
    THCONSLD(:) = 3500.0
  END IF
  IF (RDINQR('THCONWAT')) THEN
    CALL RDFREA('THCONWAT', THCONWAT, NMAT, NMAT)
  ELSE
    THCONWAT(:) = 492.0
  END IF
  IF (RDINQR('THPARX')) THEN
    CALL RDFREA('THPARX', THPARX, NMAT, NMAT)
  ELSE
    THPARX(:) = 3.55
  END IF
  IF (RDINQR('THPARA')) THEN
    CALL RDFREA('THPARA', THPARA, NMAT, NMAT)
  ELSE
    THPARA(:) = 0.75
  END IF
  IF (RDINQR('THPARB')) THEN
    CALL RDFREA('THPARB', THPARB, NMAT, NMAT)
  ELSE
    THPARB(:) = 1.2
  END IF

  RETURN
END SUBROUTINE SHT_BAS_RP


SUBROUTINE SHT_BAS_FM
!**************************************************************************************************
!  计算土壤热传导参数并重置 COF 和 RHS
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL    :: DELV1, DELV2
  REAL    :: THCON1, THCON2

  DO I = 1, NCEL
    J = IMAT(I)
    CALL U_THCON_CK(THCON(I), ONEW(I), ORES(J), OSAT(J), THCONSLD(J), THCONWAT(J), THPARX(J), &
                    THPARA(J), THPARB(J))
    CALL U_THCAP_DV(THCAP(I), ONEW(I), OSAT(J), SAND(J), CLAY(J))
    COF(I) = 0.0
    RHS(I) = 0.0
  END DO

  RETURN
END SUBROUTINE SHT_BAS_FM


SUBROUTINE SHT_BAS_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(BD))       DEALLOCATE(BD)
  IF (ALLOCATED(SAND))     DEALLOCATE(SAND)
  IF (ALLOCATED(CLAY))     DEALLOCATE(CLAY)
  IF (ALLOCATED(THCONSLD)) DEALLOCATE(THCONSLD)
  IF (ALLOCATED(THPARX))   DEALLOCATE(THPARX)
  IF (ALLOCATED(THPARA))   DEALLOCATE(THPARA)
  IF (ALLOCATED(THPARB))   DEALLOCATE(THPARB)
  IF (ALLOCATED(THCAP))    DEALLOCATE(THCAP)
  IF (ALLOCATED(THCON))    DEALLOCATE(THCON)
  IF (ALLOCATED(COF))      DEALLOCATE(COF)
  IF (ALLOCATED(RHS))      DEALLOCATE(RHS)

  RETURN
END SUBROUTINE SHT_BAS_DA

