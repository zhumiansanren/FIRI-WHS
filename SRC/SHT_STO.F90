!==================================================================================================
! 子程序包: SHT_STO (Storage Package for Soil Heat Transport Process)
! 主要功能: 模拟土壤剖面热储量的变化
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE SHT_STO_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SHT_STO_DF


SUBROUTINE SHT_STO_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(TINI)) ALLOCATE(TINI(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(TOLD)) ALLOCATE(TOLD(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(TNEW)) ALLOCATE(TNEW(NCEL), SOURCE = 0.0)

  RETURN
END SUBROUTINE SHT_STO_AL


SUBROUTINE SHT_STO_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('TINI')) THEN
    CALL RDFREA('TINI', TINI, NCEL, NCEL)
  ELSE
    TINI(:) = 20.0
  END IF

  TNEW(:) = TINI(:)
  TOLD(:) = TINI(:)

  RETURN
END SUBROUTINE SHT_STO_RP


SUBROUTINE SHT_STO_AD
!**************************************************************************************************
! 保存前一时段结束时的地温
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  TOLD(:) = TNEW(:)

  RETURN
END SUBROUTINE SHT_STO_AD


SUBROUTINE SHT_STO_FM
!**************************************************************************************************
! 计算表示土壤热储量的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I

  DO I = 1, NCEL
    COF(I) = COF(I) - DELV(I) * THCAP(I)
    RHS(I) = RHS(I) - DELV(I) * THCAP(I) * TOLD(I)
  END DO

  RETURN
END SUBROUTINE SHT_STO_FM


SUBROUTINE SHT_STO_BD
!**************************************************************************************************
! 计算表示土壤热储量的均衡项
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SHT_STO_BD


SUBROUTINE SHT_STO_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(TINI)) DEALLOCATE(TINI)
  IF (ALLOCATED(TOLD)) DEALLOCATE(TOLD)
  IF (ALLOCATED(TNEW)) DEALLOCATE(TNEW)

  RETURN
END SUBROUTINE SHT_STO_DA

