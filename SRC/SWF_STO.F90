!==================================================================================================
! 子程序包: SWF_STO (Storage Package for Soil Water Flow Process)
! 主要功能: 模拟土壤剖面水分储量的变化
! 创建日期: 2020年8月30日
! 修改日志:
!==================================================================================================


SUBROUTINE SWF_STO_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('ISWFINI')) THEN
    CALL RDSINT('ISWFINI', ISWFINI)
  ELSE
    ISWFINI = 0
  END IF

  RETURN
END SUBROUTINE SWF_STO_DF


SUBROUTINE SWF_STO_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(HINI)) ALLOCATE(HINI(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HOLD)) ALLOCATE(HOLD(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HNEW)) ALLOCATE(HNEW(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OINI)) ALLOCATE(OINI(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OOLD)) ALLOCATE(OOLD(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(ONEW)) ALLOCATE(ONEW(NCEL), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_STO_AL


SUBROUTINE SWF_STO_RP
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
  IF (RDINQR('HINI')) THEN
    CALL RDFREA('HINI', HINI, NCEL, NCEL)
  ELSE
    HINI(:) = -330.0
  END IF
  IF (RDINQR('OINI')) THEN
    CALL RDFREA('OINI', OINI, NCEL, NCEL)
  ELSE
    OINI(:) = 0.2
  END IF

  ! 根据初始条件类型换算水势和含水率
  IF (ISWFINI == 0) THEN
    DO I = 1, NCEL
      J = IMAT(I)
      CALL U_O_MVG(OINI(I), HINI(I), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    END DO
  ELSE
    DO I = 1, NCEL
      J = IMAT(I)
      CALL U_H_MVG(HINI(I), OINI(I), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    END DO
  END IF

  HOLD(:) = HINI(:)
  HNEW(:) = HINI(:)
  OOLD(:) = OINI(:)
  ONEW(:) = OINI(:)

  RETURN
END SUBROUTINE SWF_STO_RP


SUBROUTINE SWF_STO_AD
!**************************************************************************************************
! 保存前一时段结束时的水势和含水率
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  HOLD(:) = HNEW(:)
  OOLD(:) = ONEW(:)

  RETURN
END SUBROUTINE SWF_STO_AD


SUBROUTINE SWF_STO_FM
!**************************************************************************************************
! 计算表示土壤水分储量的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I

  DO I = 1, NCEL
    COF(I) = COF(I) - DELV(I) * HYCAP(I)
    RHS(I) = RHS(I) + DELV(I) * (ONEW(I) - OOLD(I) - HYCAP(I) * HNEW(I))
  END DO

  RETURN
END SUBROUTINE SWF_STO_FM


SUBROUTINE SWF_STO_BD
!**************************************************************************************************
! 计算表示土壤水分储量的均衡项
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SWF_STO_BD


SUBROUTINE SWF_STO_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(HINI)) DEALLOCATE(HINI)
  IF (ALLOCATED(HOLD)) DEALLOCATE(HOLD)
  IF (ALLOCATED(HNEW)) DEALLOCATE(HNEW)
  IF (ALLOCATED(OINI)) DEALLOCATE(OINI)
  IF (ALLOCATED(OOLD)) DEALLOCATE(OOLD)
  IF (ALLOCATED(ONEW)) DEALLOCATE(ONEW)

  RETURN
END SUBROUTINE SWF_STO_DA

