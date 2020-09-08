!==================================================================================================
! 子程序包: SST_STO (Storage Package for Soil Salt Transport Process)
! 主要功能: 模拟土壤单元体盐分储量的变化
! 创建日期: 2020年7月4日
! 修改日志:
!==================================================================================================


SUBROUTINE SST_STO_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('LADS')) THEN
    CALL RDSLOG('LADS', LADS)
  ELSE
    LADS = .FALSE.
  END IF

  RETURN
END SUBROUTINE SST_STO_DF


SUBROUTINE SST_STO_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(CADSMAX)) ALLOCATE(CADSMAX(NMAT), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(PADS))    ALLOCATE(PADS(NMAT),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(CINI))    ALLOCATE(CINI(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(CNEW))    ALLOCATE(CNEW(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(COLD))    ALLOCATE(COLD(NCEL),    SOURCE = 0.0)

  RETURN
END SUBROUTINE SST_STO_AL


SUBROUTINE SST_STO_RP
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
  IF (RDINQR('CADSMAX')) THEN
    CALL RDFREA('CADSMAX',CADSMAX, NMAT, NMAT)
  ELSE
    CADSMAX(:) = 1000.0
  END IF
  IF (RDINQR('PADS')) THEN
    CALL RDFREA('PADS', PADS, NMAT, NMAT)
  ELSE
    PADS(:) = 1.0
  END IF
  IF (RDINQR('CINI')) THEN
    CALL RDFREA('CINI', CINI, NCEL, NCEL)
  ELSE
    CINI(:) = 1.0
  END IF

  CNEW(:) = CINI(:)
  COLD(:) = CINI(:)

  RETURN
END SUBROUTINE SST_STO_RP



SUBROUTINE SST_STO_AD
!**************************************************************************************************
! 保存前一时段结束时的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  COLD(:) = CNEW(:)

  RETURN
END SUBROUTINE SST_STO_AD


SUBROUTINE SST_STO_FM
!**************************************************************************************************
! 计算表示土壤盐分储量的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL :: R1, R2             ! 吸附延迟因子

  DO I = 1, NCEL
    J = IMAT(I)
    R1 = 1.0
    R2 = 1.0
   IF (LADS) THEN
     R1 = 1.0 + BD(J) / ONEW(I) * (PADS(J) * CADSMAX(J) / (1.0 + PADS(J) * COLD(I)) ** 2.0)
     R2 = 1.0 + BD(J) / OOLD(I) * (PADS(J) * CADSMAX(J) / (1.0 + PADS(J) * COLD(I)) ** 2.0)
   END IF
    COF(I) = COF(I) - DELV(I) * ONEW(I) * R1
    RHS(I) = RHS(I) - DELV(I) * OOLD(I) * COLD(I) * R2
  END DO

  RETURN
END SUBROUTINE SST_STO_FM


SUBROUTINE SST_STO_BD
!**************************************************************************************************
! 计算表示土壤盐分储量的均衡项
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SST_STO_BD


SUBROUTINE SST_STO_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(CADSMAX)) DEALLOCATE(CADSMAX)
  IF (ALLOCATED(PADS))    DEALLOCATE(PADS)
  IF (ALLOCATED(CINI))    DEALLOCATE(CINI)
  IF (ALLOCATED(COLD))    DEALLOCATE(COLD)
  IF (ALLOCATED(CNEW))    DEALLOCATE(CNEW)

  RETURN
END SUBROUTINE SST_STO_DA

