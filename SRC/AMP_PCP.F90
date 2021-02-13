!==================================================================================================
! 子程序包: AMP_PCP (Precipitation Package for Atmosphere Process)
! 主要功能: 读取降雨数据
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE AMP_PCP_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE AMP_PCP_DF


SUBROUTINE AMP_PCP_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(WRAIN))   ALLOCATE(WRAIN(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WCANO))   ALLOCATE(WCANO(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSURF))   ALLOCATE(WSURF(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(CRAIN))   ALLOCATE(CRAIN(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSACANO)) ALLOCATE(WSACANO(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSASURF)) ALLOCATE(WSASURF(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE AMP_PCP_AL


SUBROUTINE AMP_PCP_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_MET 文件读取数据
  IF (RDINQR('WRAIN')) THEN
    CALL RDFREA('WRAIN', WRAIN, NPER, NPER)
  ELSE
    WRAIN(:) = 0.0
  END IF
  IF (RDINQR('CRAIN')) THEN
    CALL RDFREA('CRAIN', CRAIN, NPER, NPER)
  ELSE
    CRAIN(:) = 0.0
  END IF

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('WRAIN')) THEN
    CALL RDFREA('WRAIN', WRAIN, NPER, NPER)
  END IF

  IF (RDINQR('CRAIN')) THEN
    CALL RDFREA('CRAIN', CRAIN, NPER, NPER)
  END IF

  RETURN
END SUBROUTINE AMP_PCP_RP


SUBROUTINE AMP_PCP_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (LCROP) THEN

    ! 将降雨量添加到 WCANO
    WCANO(JPER) = WCANO(JPER) + WRAIN(JPER)

    ! 将降雨中的盐量添加到 WSACANO
    IF (LSALT) WSACANO(JPER) = WSACANO(JPER) + WRAIN(JPER) * CRAIN(JPER)
  ELSE
    ! 将降雨量添加到 WSURF
    WSURF(JPER) = WSURF(JPER) + WRAIN(JPER)

    ! 将降雨中的盐量添加到 WSASURF
    IF (LSALT) WSASURF(JPER) = WSASURF(JPER) + WRAIN(JPER) * CRAIN(JPER)

  END IF

  RETURN
END SUBROUTINE AMP_PCP_ST


SUBROUTINE AMP_PCP_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(WRAIN))   DEALLOCATE(WRAIN)
  IF (ALLOCATED(WCANO))   DEALLOCATE(WCANO)
  IF (ALLOCATED(WSURF))   DEALLOCATE(WSURF)
  IF (ALLOCATED(CRAIN))   DEALLOCATE(CRAIN)
  IF (ALLOCATED(WSACANO)) DEALLOCATE(WSACANO)
  IF (ALLOCATED(WSASURF)) DEALLOCATE(WSASURF)

  RETURN
END SUBROUTINE AMP_PCP_DA

