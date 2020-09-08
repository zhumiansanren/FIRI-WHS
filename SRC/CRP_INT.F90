!==================================================================================================
! 子程序包: CRP_INT (Interception Package for Soil Salt Transport Process)
! 主要功能: 模拟冠层截留作用, 计算冠层截留的水分和盐分以及贯穿冠层达到地表的水分和盐分
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE CRP_INT_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('PINT')) THEN
    CALL RDSREA('PINT', PINT)
  ELSE
    PINT = 0.025
  END IF

  RETURN
END SUBROUTINE CRP_INT_DF


SUBROUTINE CRP_INT_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(LAI))     ALLOCATE(LAI(NPER),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(SCF))     ALLOCATE(SCF(NPER),     SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WCANO))   ALLOCATE(WCANO(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WINT))    ALLOCATE(WINT(NPER),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSURF))   ALLOCATE(WSURF(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSACANO)) ALLOCATE(WSACANO(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSASURF)) ALLOCATE(WSASURF(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE CRP_INT_AL


SUBROUTINE CRP_INT_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE CRP_INT_RP


SUBROUTINE CRP_INT_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
 
  ! 计算冠层截留水量
  WINT(JPER) = 0.0
  IF ((LAI(JPER) > 0.0) .AND. (WCANO(JPER) > 0.0)) THEN
    WINT(JPER) = PINT * LAI(JPER) * (1.0 - 1.0 / (1.0 + SCF(JPER) * WCANO(JPER) / PINT / LAI(JPER)))
  END IF
  
  ! 计算贯穿冠层达到地表的水量
  WSURF(JPER) = WSURF(JPER) + WCANO(JPER) - WINT(JPER)

  ! 计算贯穿冠层达到地表的盐分量
  IF ((LSALT) .AND. (WCANO(JPER) > 0.0)) THEN
    WSASURF(JPER) = WSASURF(JPER) + WSACANO(JPER) - WINT(JPER) * WSACANO(JPER) / WCANO(JPER)
  END IF

  RETURN
END SUBROUTINE CRP_INT_ST


SUBROUTINE CRP_INT_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(LAI))     DEALLOCATE(LAI)
  IF (ALLOCATED(SCF))     DEALLOCATE(SCF)
  IF (ALLOCATED(WCANO))   DEALLOCATE(WCANO)
  IF (ALLOCATED(WINT))    DEALLOCATE(WINT)
  IF (ALLOCATED(WSURF))   DEALLOCATE(WSURF)
  IF (ALLOCATED(WSACANO)) DEALLOCATE(WSACANO)
  IF (ALLOCATED(WSASURF)) DEALLOCATE(WSASURF)

  RETURN
END SUBROUTINE CRP_INT_DA

