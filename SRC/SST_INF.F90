!==================================================================================================
! 子程序包: SST_INF (Infiltration Package for Soil Salt Transport Process)
! 主要功能: 模拟地表盐分入渗
! 创建日期: 2020年7月24日
! 修改日志: 
!==================================================================================================


SUBROUTINE SST_INF_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE
  
  RETURN
END SUBROUTINE SST_INF_DF


SUBROUTINE SST_INF_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
  
  IF (.NOT. ALLOCATED(CSURF))   ALLOCATE(CSURF(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSASURF)) ALLOCATE(WSASURF(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSARNOF)) ALLOCATE(WSARNOF(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSAINF))  ALLOCATE(WSAINF(NPER),  SOURCE = 0.0)

  RETURN
END SUBROUTINE SST_INF_AL


SUBROUTINE SST_INF_FM
!**************************************************************************************************
! 计算表示地表盐分入渗的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
  
  IF (QINF > 0.0) THEN
    COF(1) = COF(1) + 0.5 * QINF
    RHS(1) = RHS(1) - 0.5 * WSASURF(JPER) * DELT
  END IF
  
  RETURN
END SUBROUTINE SST_INF_FM


SUBROUTINE SST_INF_BD
!**************************************************************************************************
! 计算表示地表盐分入渗的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
  
  WSARNOF(JPER) = WSARNOF(JPER) + QSARNOF
  WSAINF(JPER)  = WSAINF(JPER)  + QSAINF

  RETURN
END SUBROUTINE SST_INF_BD


SUBROUTINE SST_INF_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
  
  IF (ALLOCATED(CSURF))   DEALLOCATE(CSURF)
  IF (ALLOCATED(WSASURF)) DEALLOCATE(WSASURF)
  IF (ALLOCATED(WSARNOF)) DEALLOCATE(WSARNOF)
  IF (ALLOCATED(WSAINF))  DEALLOCATE(WSAINF)
  
  RETURN
END SUBROUTINE SST_INF_DA
