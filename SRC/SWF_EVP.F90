!==================================================================================================
! 子程序包: SWF_EVP (Evaporation Pacakge for Soil Water Flow Process)
! 主要功能: 模拟地表蒸发
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE SWF_EVP_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('HATM')) THEN
    CALL RDSREA('HATM', HATM)
  ELSE
    HATM = -10000.0
  END IF

  RETURN
END SUBROUTINE SWF_EVP_DF


SUBROUTINE SWF_EVP_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(SCF))    ALLOCATE(SCF(NPER),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WEPOT))  ALLOCATE(WEPOT(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WEACT))  ALLOCATE(WEACT(NPER),  SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_EVP_AL


SUBROUTINE SWF_EVP_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_MET 文件读取数据
  CALL RDINIT(U_IN_MET, U_LOG, S_IN_MET)
  IF (RDINQR('WEPOT')) THEN
    CALL RDFREA('WEPOT', WEPOT, NPER, NPER)
  END IF

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('WEPOT')) THEN
    CALL RDFREA('WEPOT', WEPOT, NPER, NPER)
  END IF

  RETURN
END SUBROUTINE SWF_EVP_RP


SUBROUTINE SWF_EVP_FM
!**************************************************************************************************
! 计算表示地表蒸发的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL    :: HYCONATM, HYCON1
  REAL    :: DELV1
  REAL    :: QEMAX                ! 当前时段的地表最大蒸发量 (cm)

  I = 1
  J = IMAT(I)

  QEACT = 0.0
  CALL U_HYCON_MVG(HYCONATM, HATM, HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
  HYCON1 = (HYCONATM + HYCON(I)) * 0.5
  DELV1 = DELV(I) * 0.5
  QEMAX = DELT * HYCON1 * ((HNEW(I) - HATM) / DELV1 - 1.0)
  QEACT = MIN(DELT * WEPOT(JPER), QEMAX)
  IF (QEACT < 0.0) QEACT = 0.0
  IF (HPONDNEW > 0.0) THEN ! 如果地表积水, 地表蒸发为水面蒸发
    HPONDNEW = HPONDNEW - QEACT
  ELSE ! 否则, 地表蒸发为土壤蒸发
    RHS(I) = RHS(I) + QEACT
  END IF

  RETURN
END SUBROUTINE SWF_EVP_FM


SUBROUTINE SWF_EVP_BD
!**************************************************************************************************
! 计算表示地表蒸发的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  WEACT(JPER) = WEACT(JPER) + QEACT

  RETURN
END SUBROUTINE SWF_EVP_BD


SUBROUTINE SWF_EVP_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(SCF))     DEALLOCATE(SCF)
  IF (ALLOCATED(WEPOT))   DEALLOCATE(WEPOT)
  IF (ALLOCATED(WEACT))   DEALLOCATE(WEACT)

  RETURN
END SUBROUTINE SWF_EVP_DA

