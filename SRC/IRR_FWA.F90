!==================================================================================================
! 子程序包: IRR_FWA (Fixed Water Application Package for Irrigation Process)
! 主要功能: 模拟固定灌溉制度
! 创建日期: 2020年7月13日
! 修改日志:
!==================================================================================================


SUBROUTINE IRR_FWA_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('IIRRCAT')) THEN
    CALL RDSINR('IIRRCAT', 0, 2, IIRRCAT)
  ELSE
    IIRRCAT = 1
  END IF
  IF (RDINQR('LIRRINT')) THEN
    CALL RDSLOG('LIRRINT', LIRRINT)
  ELSE
    LIRRINT = .FALSE.
  END IF

  RETURN
END SUBROUTINE IRR_FWA_DF


SUBROUTINE IRR_FWA_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
  
  IF (.NOT. ALLOCATED(WIRR))    ALLOCATE(WIRR(NPER),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WCANO))   ALLOCATE(WCANO(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSURF))   ALLOCATE(WSURF(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(CIRR))    ALLOCATE(CIRR(NPER),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSACANO)) ALLOCATE(WSACANO(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSASURF)) ALLOCATE(WSASURF(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE IRR_FWA_AL


SUBROUTINE IRR_FWA_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR
  
  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('WIRR')) THEN
    CALL RDFREA('WIRR', WIRR, NPER, NPER)
  ELSE
    WIRR(:) = 0.0
  END IF
  IF (RDINQR('CIRR')) THEN
    CALL RDFREA('CIRR', CIRR, NPER, NPER)
  ELSE
    CIRR(:) = 0.0
  END IF

  RETURN
END SUBROUTINE IRR_FWA_RP


SUBROUTINE IRR_FWA_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (IIRRCAT == 1) THEN
    IF (LCROP .AND. LIRRINT) THEN

      ! 将灌水量添加到 WCANO
      WCANO(JPER) = WCANO(JPER) + WIRR(JPER)

      ! 将灌溉水中的盐量添加到 WSACANO
      IF (LSALT) WSACANO(JPER) = WSACANO(JPER) + WIRR(JPER) * CIRR(JPER)

    ELSE
      ! 将灌水量添加到 WSURF
      WSURF(JPER) = WSURF(JPER) + WIRR(JPER)

      ! 将灌溉水中的盐量添加到 WSASURF
      IF (LSALT) WSASURF(JPER) = WSASURF(JPER) + WIRR(JPER) * CIRR(JPER)

    END IF
  END IF

  RETURN
END SUBROUTINE IRR_FWA_ST


SUBROUTINE IRR_FWA_FM
!**************************************************************************************************
! 计算方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (IIRRCAT == 2) THEN
    IF ((ICELIRR > 0) .AND. (ICELIRR <= NCEL)) &
      RHS(ICELIRR) = RHS(ICELIRR) + DELT * WIRR(JPER)
  END IF

  RETURN
END SUBROUTINE IRR_FWA_FM


SUBROUTINE IRR_FWA_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(WIRR))    DEALLOCATE(WIRR)
  IF (ALLOCATED(WCANO))   DEALLOCATE(WCANO)
  IF (ALLOCATED(WSURF))   DEALLOCATE(WSURF)
  IF (ALLOCATED(CIRR))    DEALLOCATE(CIRR)
  IF (ALLOCATED(WSACANO)) DEALLOCATE(WSACANO)
  IF (ALLOCATED(WSASURF)) DEALLOCATE(WSASURF)

  RETURN
END SUBROUTINE IRR_FWA_DA

