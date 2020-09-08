!==================================================================================================
! 子程序包: IRR_BAS (Basic Package for Irrigation Process)
! 主要功能: 为模型提供灌溉数据
! 创建日期: 2020年7月13日
! 修改日志:
!==================================================================================================


SUBROUTINE IRR_BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('IIRRSTG')) THEN
    CALL RDSINT('IIRRSTG', IIRRSTG)
  ELSE
    IIRRSTG = 1
  END IF

  L_IRR_FWA = (IIRRSTG == 1)
  L_IRR_OWA = (IIRRSTG == 2)

  RETURN
END SUBROUTINE IRR_BAS_DF


SUBROUTINE IRR_BAS_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
  
  IF (.NOT. ALLOCATED(WIRR)) ALLOCATE(WIRR(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(CIRR)) ALLOCATE(CIRR(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE IRR_BAS_AL


SUBROUTINE IRR_BAS_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE IRR_BAS_RP


SUBROUTINE IRR_BAS_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE IRR_BAS_ST


SUBROUTINE IRR_BAS_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(WIRR)) DEALLOCATE(WIRR)
  IF (ALLOCATED(CIRR)) DEALLOCATE(CIRR)

  RETURN
END SUBROUTINE IRR_BAS_DA

