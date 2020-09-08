!==================================================================================================
! 子程序包: CRP_BAS (Basic Package for Crop Process)
! 主要功能: 为模型提供覆盖度, 叶面积指数, 株高, 根系深度等数据
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE CRP_BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('ISCF')) THEN
    CALL RDSINT('ISCF', ISCF)
  ELSE
    ISCF = 1
  END IF
  IF (RDINQR('EXTINC')) THEN
    CALL RDSREA('EXTINC', EXTINC)
  ELSE
    EXTINC = 0.39
  END IF
  IF (RDINQR('ALBEDO')) THEN
    CALL RDSREA('ALBEDO', ALBEDO)
  ELSE
    ALBEDO = 0.23
  END IF

  LRWU = .TRUE.

  L_CRP_PET = .TRUE.
  L_CRP_INT = .TRUE.

  RETURN
END SUBROUTINE CRP_BAS_DF


SUBROUTINE CRP_BAS_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE
  
  IF (.NOT. ALLOCATED(SCF))    ALLOCATE(SCF(NPER),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(LAI))    ALLOCATE(LAI(NPER),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RDEPTH)) ALLOCATE(RDEPTH(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(CRPHGT)) ALLOCATE(CRPHGT(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE CRP_BAS_AL


SUBROUTINE CRP_BAS_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_CRP 文件读取数据
  CALL RDINIT(U_IN_CRP, U_LOG, S_IN_CRP)
  IF (RDINQR('LAI')) THEN
    CALL RDFREA('LAI', LAI, NPER, NPER)
  ELSE
    LAI(:) = 2.5
  END IF
  IF (RDINQR('SCF')) THEN
    CALL RDFREA('SCF', SCF, NPER, NPER)
  ELSE
    SCF(:) = 0.5
  END IF
  IF (RDINQR('RDEPTH')) THEN
    CALL RDFREA('RDEPTH', RDEPTH, NPER, NPER)
  ELSE
    RDEPTH(:) = 0.5 * LOCV(NCEL, 2)
  END IF
  IF (RDINQR('CRPHGT')) THEN
    CALL RDFREA('CRPHGT', CRPHGT, NPER, NPER)
  ELSE
    CRPHGT(:) = 50.0
  END IF

  RETURN
END SUBROUTINE CRP_BAS_RP


SUBROUTINE CRP_BAS_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 换算 LAI 和 SCF
  IF (ISCF == 0) THEN
    LAI(JPER) = - LOG(1.0 - SCF(JPER)) / EXTINC
  ELSE
    SCF(JPER) = 1.0 - EXP(- EXTINC * LAI(JPER))
  END IF

  RETURN
END SUBROUTINE CRP_BAS_ST


SUBROUTINE CRP_BAS_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(SCF)) DEALLOCATE(SCF)
  IF (ALLOCATED(LAI)) DEALLOCATE(LAI)
  IF (ALLOCATED(CRPHGT)) DEALLOCATE(CRPHGT)
  IF (ALLOCATED(RDEPTH)) DEALLOCATE(RDEPTH)

  RETURN
END SUBROUTINE CRP_BAS_DA

