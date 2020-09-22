!==================================================================================================
! 子程序包: BAS (Basic Package)
! 主要功能: 读取模型的基本设置, 计算时段步长, 输出模拟结果
! 创建日期: 2020年7月3日
! 修改日志:
!
!==================================================================================================


SUBROUTINE BAS_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 打开 S_LOG 文件, 输出标识信息
  OPEN(UNIT = U_LOG, FILE = S_LOG)
  WRITE(U_LOG, "(A)") APPNAM
  WRITE(U_LOG, "()")

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('PRJNAM')) THEN
    CALL RDSCHA('PRJNAM', PRJNAM)
  ELSE
    PRJNAM = 'This is a default project.'
  END IF
  IF (RDINQR('IRUNMOD')) THEN
    CALL RDSINT('IRUNMOD', IRUNMOD)
  ELSE
    IRUNMOD = 1
  END IF
  IF (RDINQR('NMAT')) THEN
    CALL RDSINT('NMAT', NMAT)
  ELSE
    NMAT = 1
  END IF
  IF (RDINQR('NCEL')) THEN
    CALL RDSINT('NCEL', NCEL)
  ELSE
    NCEL = 100
  END IF
  IF (RDINQR('NPER')) THEN
    CALL RDSINT('NPER', NPER)
  ELSE
    NPER = 10
  END IF
  IF (RDINQR('LMETEO')) THEN
    CALL RDSLOG('LMETEO', LMETEO)
  ELSE
    LMETEO = .FALSE.
  END IF
  IF (RDINQR('LCROP')) THEN
    CALL RDSLOG('LCROP', LCROP)
  ELSE
    LCROP = .FALSE.
  END IF
  IF (RDINQR('LIRRI')) THEN
    CALL RDSLOG('LIRRI', LIRRI)
  ELSE
    LIRRI = .FALSE.
  END IF
  IF (RDINQR('LHEAT')) THEN
    CALL RDSLOG('LHEAT', LHEAT)
  ELSE
    LHEAT = .FALSE.
  END IF
  IF (RDINQR('LSALT')) THEN
    CALL RDSLOG('LSALT', LSALT)
  ELSE
    LSALT = .FALSE.
  END IF
  IF (RDINQR('LRWU')) THEN
    CALL RDSLOG('LRWU', LRWU)
  ELSE
    LRWU = .FALSE.
  END IF
   IF (RDINQR('LSSK')) THEN
    CALL RDSLOG('LSSK', LSSK)
  ELSE
    LSSK = .FALSE.
  END IF
  IF (RDINQR('S_IN_SOL')) THEN
    CALL RDSCHA('S_IN_SOL', S_IN_SOL)
  ELSE
    S_IN_SOL = S_IN_BAS
  END IF
  IF (RDINQR('S_IN_PER')) THEN
    CALL RDSCHA('S_IN_PER', S_IN_PER)
  ELSE
    S_IN_PER = S_IN_BAS
  END IF
  IF (RDINQR('S_IN_MET')) THEN
    CALL RDSCHA('S_IN_MET', S_IN_MET)
  ELSE
    S_IN_MET = S_IN_BAS
  END IF
  IF (RDINQR('S_IN_CRP')) THEN
    CALL RDSCHA('S_IN_CRP', S_IN_CRP)
  ELSE
    S_IN_CRP = S_IN_BAS
  END IF
  IF (RDINQR('S_OUT_LST')) THEN
    CALL RDSCHA('S_OUT_LST', S_OUT_LST)
  ELSE
    S_OUT_LST = 'LIST.OUT'
  END IF
  IF (RDINQR('S_OUT_PER')) THEN
    CALL RDSCHA('S_OUT_PER', S_OUT_PER)
  ELSE
    S_OUT_PER = 'PERIOD.OUT'
  END IF
  IF (RDINQR('S_OUT_OBS')) THEN
    CALL RDSCHA('S_OUT_OBS', S_OUT_OBS)
  ELSE
    S_OUT_OBS = 'OBS.OUT'
  END IF

  L_SWF_BAS = .TRUE.
  L_AMP_BAS = LMETEO
  L_CRP_BAS = LCROP
  L_IRR_BAS = LIRRI
  L_SHT_BAS = LHEAT
  L_SST_BAS = LSALT

  RETURN
END SUBROUTINE BAS_DF


SUBROUTINE BAS_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(PERLEN)) ALLOCATE(PERLEN(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(NSTP))   ALLOCATE(NSTP(NPER),    SOURCE = 0)
  IF (.NOT. ALLOCATED(TSMLT))  ALLOCATE(TSMLT(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(IPRNT))  ALLOCATE(IPRNT(NPER),   SOURCE = 0)
  IF (.NOT. ALLOCATED(IMAT))   ALLOCATE(IMAT(NCEL),    SOURCE = 0)
  IF (.NOT. ALLOCATED(IOBS))   ALLOCATE(IOBS(NCEL),    SOURCE = 0)
  IF (.NOT. ALLOCATED(DELV))   ALLOCATE(DELV(NCEL),    SOURCE = 0.0)
  IF (.NOT. ALLOCATED(LOCV))   ALLOCATE(LOCV(NCEL, 3), SOURCE = 0.0)

  RETURN
END SUBROUTINE BAS_AL


SUBROUTINE BAS_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 局部变量
  INTEGER :: I

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('PERLEN')) THEN
    CALL RDFREA('PERLEN', PERLEN, NPER, NPER)
  ELSE
    PERLEN(:) = 1.0
  END IF
  IF (RDINQR('NSTP')) THEN
    CALL RDFINT('NSTP', NSTP, NPER, NPER)
  ELSE
    NSTP(:) = 30
  END IF
  IF (RDINQR('TSMLT')) THEN
    CALL RDFREA('TSMLT', TSMLT, NPER, NPER)
  ELSE
    TSMLT(:) = 1.05
  END IF
  IF (RDINQR('IPRNT')) THEN
    CALL RDFINT('IPRNT', IPRNT, NPER, NPER)
  ELSE
    IPRNT(:) = 1
  END IF

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('IMAT')) THEN
    CALL RDFINT('IMAT', IMAT, NCEL, NCEL)
  ELSE
    IMAT = 1
  END IF
  IF (RDINQR('DELV')) THEN
    CALL RDFREA('DELV', DELV, NCEL, NCEL)
  ELSE
    DELV(:) = 1.0
  END IF
  IF (RDINQR('IOBS')) THEN
    CALL RDFINT('IOBS', IOBS, NCEL, NCEL)
  ELSE
    IOBS(:) = 0
  END IF

  ! 打开 S_OUT_LST 文件, 输出标识信息
  OPEN(UNIT = U_OUT_LST, FILE = S_OUT_LST)
  WRITE(U_OUT_LST, "(A)") APPNAM
  WRITE(U_OUT_LST, "()")
  WRITE(U_OUT_LST, "('项目名称: ', A)") TRIM(PRJNAM)
  WRITE(U_OUT_LST, "('应力期数量: ', I12)") NPER
  WRITE(U_OUT_LST, "('单元体数量: ', I12)") NCEL
  WRITE(U_OUT_LST, "()")

  ! 打开 S_OUT_PER 文件, 输出标识信息
  OPEN(UNIT = U_OUT_PER, FILE = S_OUT_PER)
  WRITE(U_OUT_PER, "(A)") APPNAM
  WRITE(U_OUT_PER, "()")
  WRITE(U_OUT_PER, "('项目名称: ', A)") TRIM(PRJNAM)
  WRITE(U_OUT_PER, "('应力期数量: ', I12)") NPER
  WRITE(U_OUT_PER, "()")

  ! 打开 S_OUT_OBS 文件, 输出标识信息
  OPEN(UNIT = U_OUT_OBS, FILE = S_OUT_OBS)
  WRITE(U_OUT_OBS, '(A)') APPNAM
  WRITE(U_OUT_OBS, "()")
  WRITE(U_OUT_OBS, "('项目名称: ', A)") TRIM(PRJNAM)
  WRITE(U_OUT_OBS, "('应力期数量: ', I12)") NPER
  WRITE(U_OUT_OBS, "()")

  ! 在屏幕上显示标识信息
  IF (IRUNMOD == 1) THEN
    WRITE(*, "(A)") APPNAM
    WRITE(*, "()")
    WRITE(*, "('项目名称: ', A)") TRIM(PRJNAM)
    WRITE(*, "('应力期数量: ', I12)") NPER
    WRITE(*, "()")
  END IF

  ! 计算土壤单元体顶部, 中心和底部的位置
  LOCV(1, 1) = 0.0
  LOCV(1, 2) = DELV(1) * 0.5
  LOCV(1, 3) = DELV(1)
  DO I = 2, NCEL
    LOCV(I, 1) = LOCV(I - 1, 1) + DELV(I)
    LOCV(I, 2) = LOCV(I - 1, 2) + DELV(I)
    LOCV(I, 3) = LOCV(I - 1, 3) + DELV(I)
  END DO

  ! 初始化模拟总时长
  TOTIM = 0.0

  ! 初始化错误代码
  IERRCOD = 0

  RETURN
END SUBROUTINE BAS_RP


SUBROUTINE BAS_ST
!**************************************************************************************************
! 计算当前应力期的初始时段步长
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 计算当前应力期的初始时段步长
  DELT = PERLEN(JPER) / FLOAT(NSTP(JPER))
  IF (TSMLT(JPER) /= 1.0) THEN
    DELT = PERLEN(JPER) * (1.0 - TSMLT(JPER)) / (1.0 - TSMLT(JPER) ** FLOAT(NSTP(JPER)))
  END IF

  ! 将当前应力期的序号输出到日志文件并在屏幕上显示
  WRITE(U_LOG, "('  应力期: ', I12)") JPER
  IF (IRUNMOD == 1) WRITE(*, "('  应力期: ', I12)") JPER

  RETURN
END SUBROUTINE BAS_ST


SUBROUTINE BAS_AD
!**************************************************************************************************
! 计算当前时段步长和模拟总时长
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 计算当前时段步长
  IF (JSTP /= 1) DELT = DELT * TSMLT(JPER)

  ! 计算模拟总时长
  TOTIM = TOTIM + DELT

  RETURN
END SUBROUTINE BAS_AD


SUBROUTINE BAS_OT
!**************************************************************************************************
! 输出模拟结果
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I

  IF (JPER == 1) THEN

    ! 将标题输出到 S_OUT_LST 文件
    WRITE(U_OUT_LST, "('      PERIOD')", ADVANCE = 'NO')
    WRITE(U_OUT_LST, "('        TIME')", ADVANCE = 'NO')
    WRITE(U_OUT_LST, "('        CELL')", ADVANCE = 'NO')
    WRITE(U_OUT_LST, "('       DEPTH')", ADVANCE = 'NO')
    IF (ALLOCATED(HNEW))    WRITE(U_OUT_LST, "('        HEAD')", ADVANCE = 'NO')
    IF (ALLOCATED(ONEW))    WRITE(U_OUT_LST, "('       MOIST')", ADVANCE = 'NO')
    IF (ALLOCATED(TNEW))    WRITE(U_OUT_LST, "('        TEMP')", ADVANCE = 'NO')
    IF (ALLOCATED(CNEW))    WRITE(U_OUT_LST, "('        CONC')", ADVANCE = 'NO')
    IF (ALLOCATED(WRWUACT)) WRITE(U_OUT_LST, "('         RWU')", ADVANCE = 'NO')
    IF (ALLOCATED(WSS))     WRITE(U_OUT_LST, "('         WSS')", ADVANCE = 'NO')
    WRITE(U_OUT_LST, "()", ADVANCE = 'YES')

    ! 将标题输出到 S_OUT_PER 文件
    WRITE(U_OUT_PER, "('      PERIOD')", ADVANCE = 'NO')
    WRITE(U_OUT_PER, "('      PERLEN')", ADVANCE = 'NO')
    WRITE(U_OUT_PER, "('        TIME')", ADVANCE = 'NO')
    IF (ALLOCATED(ONEW))   WRITE(U_OUT_PER, "('     STORAGE')", ADVANCE = 'NO')
    IF (ALLOCATED(WBOT))   WRITE(U_OUT_PER, "('        WBOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WBOT))   WRITE(U_OUT_PER, "('      C_WBOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTOP))   WRITE(U_OUT_PER, "('        WTOP')", ADVANCE = 'NO')
    IF (ALLOCATED(WTOP))   WRITE(U_OUT_PER, "('      C_WTOP')", ADVANCE = 'NO')
    IF (ALLOCATED(WIRR))   WRITE(U_OUT_PER, "('        WIRR')", ADVANCE = 'NO')
    IF (ALLOCATED(WIRR))   WRITE(U_OUT_PER, "('      C_WIRR')", ADVANCE = 'NO')
    IF (ALLOCATED(WRAIN))  WRITE(U_OUT_PER, "('       WRAIN')", ADVANCE = 'NO')
    IF (ALLOCATED(WRAIN))  WRITE(U_OUT_PER, "('     C_WRAIN')", ADVANCE = 'NO')
    IF (ALLOCATED(WCANO))  WRITE(U_OUT_PER, "('       WCANO')", ADVANCE = 'NO')
    IF (ALLOCATED(WCANO))  WRITE(U_OUT_PER, "('     C_WCANO')", ADVANCE = 'NO')
    IF (ALLOCATED(WSURF))  WRITE(U_OUT_PER, "('       WSURF')", ADVANCE = 'NO')
    IF (ALLOCATED(WSURF))  WRITE(U_OUT_PER, "('     C_WSURF')", ADVANCE = 'NO')
    IF (ALLOCATED(WINT))   WRITE(U_OUT_PER, "('        WINT')", ADVANCE = 'NO')
    IF (ALLOCATED(WINT))   WRITE(U_OUT_PER, "('      C_WINT')", ADVANCE = 'NO')
    IF (ALLOCATED(WRNOF))  WRITE(U_OUT_PER, "('       WRNOF')", ADVANCE = 'NO')
    IF (ALLOCATED(WRNOF))  WRITE(U_OUT_PER, "('     C_WRNOF')", ADVANCE = 'NO')
    IF (ALLOCATED(WINF))   WRITE(U_OUT_PER, "('        WINF')", ADVANCE = 'NO')
    IF (ALLOCATED(WINF))   WRITE(U_OUT_PER, "('      C_WINF')", ADVANCE = 'NO')
    IF (ALLOCATED(WETPOT)) WRITE(U_OUT_PER, "('      WETPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WETPOT)) WRITE(U_OUT_PER, "('    C_WETPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WEPOT))  WRITE(U_OUT_PER, "('       WEPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WEPOT))  WRITE(U_OUT_PER, "('     C_WEPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WEACT))  WRITE(U_OUT_PER, "('       WEACT')", ADVANCE = 'NO')
    IF (ALLOCATED(WEACT))  WRITE(U_OUT_PER, "('     C_WEACT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTPOT))  WRITE(U_OUT_PER, "('       WTPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTPOT))  WRITE(U_OUT_PER, "('     C_WTPOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTACT))  WRITE(U_OUT_PER, "('       WTACT')", ADVANCE = 'NO')
    IF (ALLOCATED(WTACT))  WRITE(U_OUT_PER, "('     C_WTACT')", ADVANCE = 'NO')

    IF (ALLOCATED(WSABOT)) WRITE(U_OUT_PER, "('      WSABOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WSABOT)) WRITE(U_OUT_PER, "('    C_WSABOT')", ADVANCE = 'NO')
    IF (ALLOCATED(WSATOP)) WRITE(U_OUT_PER, "('      WSATOP')", ADVANCE = 'NO')
    IF (ALLOCATED(WSATOP)) WRITE(U_OUT_PER, "('    C_WSATOP')", ADVANCE = 'NO')

    WRITE(U_OUT_PER, "()", ADVANCE = 'YES')

    ! 将标题输出到 S_OUT_OBS 文件
    WRITE(U_OUT_OBS, "('        IOBS')", ADVANCE = 'NO')
    WRITE(U_OUT_OBS, "('        CELL')", ADVANCE = 'NO')
    WRITE(U_OUT_OBS, "('       DEPTH')", ADVANCE = 'NO')
    WRITE(U_OUT_OBS, "('      PERIOD')", ADVANCE = 'NO')
    WRITE(U_OUT_OBS, "('      PERLEN')", ADVANCE = 'NO')
    WRITE(U_OUT_OBS, "('        TIME')", ADVANCE = 'NO')
    IF (ALLOCATED(HNEW)) WRITE(U_OUT_OBS, "('        HEAD')", ADVANCE = 'NO')
    IF (ALLOCATED(ONEW)) WRITE(U_OUT_OBS, "('       MOIST')", ADVANCE = 'NO')
    IF (ALLOCATED(TNEW)) WRITE(U_OUT_OBS, "('        TEMP')", ADVANCE = 'NO')
    IF (ALLOCATED(CNEW)) WRITE(U_OUT_OBS, "('        CONC')", ADVANCE = 'NO')
    WRITE(U_OUT_OBS, "()", ADVANCE = 'YES')

    ! 将初始条件数据输出到 S_OUT_LST 文件
    DO I = 1, NCEL
      WRITE(U_OUT_LST, "(I12)",   ADVANCE = 'NO') 0
      WRITE(U_OUT_LST, "(F12.4)", ADVANCE = 'NO') 0.0
      WRITE(U_OUT_LST, "(I12)",   ADVANCE = 'NO') I
      IF (ALLOCATED(LOCV))    WRITE(U_OUT_LST, "(F12.2)",  ADVANCE = 'NO') LOCV(I, 2)
      IF (ALLOCATED(HINI))    WRITE(U_OUT_LST, "(ES12.3)", ADVANCE = 'NO') HINI(I)
      IF (ALLOCATED(OINI))    WRITE(U_OUT_LST, "(F12.4)",  ADVANCE = 'NO') OINI(I)
      IF (ALLOCATED(TINI))    WRITE(U_OUT_LST, "(F12.4)",  ADVANCE = 'NO') TINI(I)
      IF (ALLOCATED(CINI))    WRITE(U_OUT_LST, "(ES12.3)", ADVANCE = 'NO') CINI(I)
      IF (ALLOCATED(WRWUACT)) WRITE(U_OUT_LST, "(ES12.3)", ADVANCE = 'NO') QRWUACT(I) / DELT
      IF (ALLOCATED(WSS))     WRITE(U_OUT_LST, "(ES12.3)", ADVANCE = 'NO') WSS(I)
      WRITE(U_OUT_LST, "()", ADVANCE = 'YES')
    END DO

  END IF

  ! 将模拟结果输出到 S_OUT_LST 文件
  IF ((IPRNT(JPER) == 1) .OR. (IPRNT(JPER) == 2) .OR. (IPRNT(JPER) == 3)) THEN
    DO I = 1, NCEL
      WRITE(U_OUT_LST, "(I12)",   ADVANCE = 'NO') JPER
      WRITE(U_OUT_LST, "(F12.4)", ADVANCE = 'NO') TOTIM
      WRITE(U_OUT_LST, "(I12)",   ADVANCE = 'NO') I
      IF (ALLOCATED(LOCV))    WRITE(U_OUT_LST, "(F12.2)",  ADVANCE = 'NO') LOCV(I, 2)
      IF (ALLOCATED(HNEW))    WRITE(U_OUT_LST, "(ES12.3)", ADVANCE = 'NO') HNEW(I)
      IF (ALLOCATED(ONEW))    WRITE(U_OUT_LST, "(F12.4)",  ADVANCE = 'NO') ONEW(I)
      IF (ALLOCATED(TNEW))    WRITE(U_OUT_LST, "(F12.4)",  ADVANCE = 'NO') TNEW(I)
      IF (ALLOCATED(CNEW))    WRITE(U_OUT_LST, "(ES12.3)", ADVANCE = 'NO') CNEW(I)
      IF (ALLOCATED(WRWUACT)) WRITE(U_OUT_LST, "(ES12.3)", ADVANCE = 'NO') QRWUACT(I) / DELT
      IF (ALLOCATED(WSS))     WRITE(U_OUT_LST, "(ES12.3)", ADVANCE = 'NO') WSS(I)
      WRITE(U_OUT_LST, "()", ADVANCE = 'YES')
    END DO
  END IF

  ! 将模拟结果输出到 S_OUT_PER 文件
  IF ((IPRNT(JPER) == 1) .OR. (IPRNT(JPER) == 2) .OR. (IPRNT(JPER) == 4)) THEN
    WRITE(U_OUT_PER, "(I12)",   ADVANCE = 'NO') JPER
    WRITE(U_OUT_PER, "(F12.4)", ADVANCE = 'NO') PERLEN(JPER)
    WRITE(U_OUT_PER, "(F12.4)", ADVANCE = 'NO') TOTIM
    IF (ALLOCATED(ONEW))   WRITE(U_OUT_PER, "(F12.4)",  ADVANCE = 'NO') SUM(ONEW(:) * DELV(:))
    IF (ALLOCATED(WBOT))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WBOT(JPER)
    IF (ALLOCATED(WBOT))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WBOT(1:JPER))
    IF (ALLOCATED(WTOP))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WTOP(JPER)
    IF (ALLOCATED(WTOP))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WTOP(1:JPER))
    IF (ALLOCATED(WIRR))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WIRR(JPER)
    IF (ALLOCATED(WIRR))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WIRR(1:JPER))
    IF (ALLOCATED(WRAIN))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WRAIN(JPER)
    IF (ALLOCATED(WRAIN))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WRAIN(1:JPER))
    IF (ALLOCATED(WCANO))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WCANO(JPER)
    IF (ALLOCATED(WCANO))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WCANO(1:JPER))
    IF (ALLOCATED(WSURF))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WSURF(JPER)
    IF (ALLOCATED(WSURF))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WSURF(1:JPER))
    IF (ALLOCATED(WINT))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WINT(JPER)
    IF (ALLOCATED(WINT))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WINT(1:JPER))
    IF (ALLOCATED(WRNOF))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WRNOF(JPER)
    IF (ALLOCATED(WRNOF))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WRNOF(1:JPER))
    IF (ALLOCATED(WINF))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WINF(JPER)
    IF (ALLOCATED(WINF))   WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WINF(1:JPER))
    IF (ALLOCATED(WETPOT)) WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WETPOT(JPER)
    IF (ALLOCATED(WETPOT)) WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WETPOT(1:JPER))
    IF (ALLOCATED(WEPOT))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WEPOT(JPER)
    IF (ALLOCATED(WEPOT))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WEPOT(1:JPER))
    IF (ALLOCATED(WEACT))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WEACT(JPER)
    IF (ALLOCATED(WEACT))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WEACT(1:JPER))
    IF (ALLOCATED(WTPOT))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WTPOT(JPER)
    IF (ALLOCATED(WTPOT))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WTPOT(1:JPER))
    IF (ALLOCATED(WTACT))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WTACT(JPER)
    IF (ALLOCATED(WTACT))  WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WTACT(1:JPER))

    IF (ALLOCATED(WSABOT)) WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WSABOT(JPER)
    IF (ALLOCATED(WSABOT)) WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WSABOT(1:JPER))
    IF (ALLOCATED(WSATOP)) WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') WSATOP(JPER)
    IF (ALLOCATED(WSATOP)) WRITE(U_OUT_PER, "(ES12.3)", ADVANCE = 'NO') SUM(WSATOP(1:JPER))

    WRITE(U_OUT_PER, "()", ADVANCE = 'YES')
  END IF


  ! 将模拟结果输出到 S_OUT_OBS 文件
  IF ((IPRNT(JPER) == 1) .OR. (IPRNT(JPER) == 3) .OR. (IPRNT(JPER) == 4)) THEN
    DO I = 1, NCEL
      IF (IOBS(I) > 0) THEN
        WRITE(U_OUT_OBS, "(I12)",   ADVANCE = 'NO') IOBS(I)
        WRITE(U_OUT_OBS, "(I12)",   ADVANCE = 'NO') I
        WRITE(U_OUT_OBS, "(F12.2)", ADVANCE = 'NO') LOCV(I, 2)
        WRITE(U_OUT_OBS, "(I12)",   ADVANCE = 'NO') JPER
        WRITE(U_OUT_OBS, "(F12.4)", ADVANCE = 'NO') PERLEN(JPER)
        WRITE(U_OUT_OBS, "(F12.4)", ADVANCE = 'NO') TOTIM
        IF (ALLOCATED(HNEW)) WRITE(U_OUT_OBS, "(ES12.3)", ADVANCE = 'NO') HNEW(I)
        IF (ALLOCATED(ONEW)) WRITE(U_OUT_OBS, "(F12.4)",  ADVANCE = 'NO') ONEW(I)
        IF (ALLOCATED(TNEW)) WRITE(U_OUT_OBS, "(F12.4)",  ADVANCE = 'NO') TNEW(I)
        IF (ALLOCATED(CNEW)) WRITE(U_OUT_OBS, "(ES12.3)", ADVANCE = 'NO') CNEW(I)
        WRITE(U_OUT_OBS, "()", ADVANCE = 'YES')
      END IF
    END DO
  END IF

  RETURN
END SUBROUTINE BAS_OT


SUBROUTINE BAS_DA
!**************************************************************************************************
! 释放数组变量内存并关闭文件
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 释放数组变量内存
  IF (ALLOCATED(PERLEN)) DEALLOCATE(PERLEN)
  IF (ALLOCATED(NSTP))   DEALLOCATE(NSTP)
  IF (ALLOCATED(TSMLT))  DEALLOCATE(TSMLT)
  IF (ALLOCATED(IPRNT))  DEALLOCATE(IPRNT)
  IF (ALLOCATED(IMAT))   DEALLOCATE(IMAT)
  IF (ALLOCATED(DELV))   DEALLOCATE(DELV)
  IF (ALLOCATED(IOBS))   DEALLOCATE(IOBS)
  IF (ALLOCATED(LOCV))   DEALLOCATE(LOCV)

  ! 删除临时文件
  CALL RDDTMP(U_IN_BAS)
  CALL RDDTMP(U_IN_SOL)
  CALL RDDTMP(U_IN_SOL)
  CALL RDDTMP(U_IN_PER)
  CALL RDDTMP(U_IN_MET)
  CALL RDDTMP(U_IN_CRP)

  ! 关闭输入文件
  CLOSE(U_IN_BAS)
  CLOSE(U_IN_SOL)
  CLOSE(U_IN_PER)
  CLOSE(U_IN_MET)
  CLOSE(U_IN_CRP)

  ! 关闭输出文件
  CLOSE(U_OUT_LST)
  CLOSE(U_OUT_PER)
  CLOSE(U_OUT_OBS)
  CLOSE(U_LOG)

  RETURN
END SUBROUTINE BAS_DA
