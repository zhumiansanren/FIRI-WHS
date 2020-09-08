!==================================================================================================
! 子程序包: SST_CCA (Cell-Centered Advection Package for Soil Salt Transport Process)
! 主要功能: 模拟土壤单元体之间盐分的对流迁移
! 创建日期: 2020年7月4日
! 修改日志:
!==================================================================================================


SUBROUTINE SST_CCA_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SST_CCA_DF


SUBROUTINE SST_CCA_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(BBV)) ALLOCATE(BBV(NCEL, 2), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(DDV)) ALLOCATE(DDV(NCEL, 2), SOURCE = 0.0)

  RETURN
END SUBROUTINE SST_CCA_AL


SUBROUTINE SST_CCA_FM
!**************************************************************************************************
! 计算表示土壤单元体之间盐分对流迁移的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
  
  INTEGER :: I

  ! 第 2 至 NCEL - 1 单元体
  DO I = 2, NCEL - 1
    BBV(I, 1) = QV(I, 1) * DELV(I - 1) / (DELV(I - 1) + DELV(I))
    BBV(I, 2) = QV(I, 2) * DELV(I + 1) / (DELV(I + 1) + DELV(I))
    DDV(I, 1) = QV(I, 1) * DELV(I) / (DELV(I - 1) + DELV(I))
    DDV(I, 2) = QV(I, 2) * DELV(I) / (DELV(I + 1) + DELV(I))
  END DO
  
  ! 第 1 单元体
  I = 1
  BBV(I, 2) = QV(I, 2) * DELV(I + 1) / (DELV(I + 1) + DELV(I))
  DDV(I, 2) = QV(I, 2) * DELV(I) / (DELV(I + 1) + DELV(I))
  
  ! 第 NCEL 单元体
  I = NCEL
  BBV(I, 1) = QV(I, 1) * DELV(I - 1) / (DELV(I - 1) + DELV(I))
  DDV(I, 1) = QV(I, 1) * DELV(I) / (DELV(I - 1) + DELV(I))

  RETURN
END SUBROUTINE SST_CCA_FM


SUBROUTINE SST_CCA_BD
!**************************************************************************************************
! 计算表示土壤单元体之间盐分对流迁移的均衡项
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SST_CCA_BD


SUBROUTINE SST_CCA_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(BBV)) DEALLOCATE(BBV)
  IF (ALLOCATED(DDV)) DEALLOCATE(DDV)

  RETURN
END SUBROUTINE SST_CCA_DA

