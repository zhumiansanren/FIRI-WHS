!==================================================================================================
! 子程序包: UTL (Utility Package)
! 主要功能: 公用函数
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE U_HYCON_MVG(K, H, KS, A, N, L)
!**************************************************************************************************
! 采用 Mualem (1976) 和 van Genuchten (1980) 公式计算土壤非饱和导水率
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: K                         ! 土壤非饱和导水率 (cm d^-1)

  ! 输入参数:
  REAL, INTENT(IN)  :: H                         ! 土壤水势 (cm)
  REAL, INTENT(IN)  :: KS                        ! 土壤饱和导水率 (cm d^-1)
  REAL, INTENT(IN)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  REAL, INTENT(IN)  :: N                         ! 土壤孔径分布指数
  REAL, INTENT(IN)  :: L                         ! 土壤孔隙连结性参数

  ! 局部变量:
  REAL              :: M                         ! M = 1 - 1 / N
  REAL              :: SE                        ! 土壤有效饱和度

  K = KS
  IF (H < 0.0) THEN
    M = 1.0 - 1.0 / N
    SE = (1.0 + (- A * H) ** N) ** (- M)
    K = KS * SE ** L * (1.0 - (1.0 - SE ** (1.0 / M)) ** M) ** 2.0
  END IF

  RETURN
END SUBROUTINE U_HYCON_MVG


SUBROUTINE U_O_MVG(O, H, OR, OS, A, N)
!**************************************************************************************************
! 采用 Mualem (1976) 和 van Genuchten (1980) 公式计算土壤含水率
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: O                         ! 土壤含水率 (cm^3 cm^-3)

  ! 输入参数:
  REAL, INTENT(IN)  :: H                         ! 土壤水势 (cm)
  REAL, INTENT(IN)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  REAL, INTENT(IN)  :: N                         ! 土壤孔径分布指数

  ! 局部变量:
  REAL              :: M                         ! M = 1 - 1 / N

  O = OS
  IF (H < 0.0) THEN
    M = 1.0 - 1.0 / N
    O = OR + (OS - OR) * (1.0 + (- A * H) ** N) ** (- M)
  END IF

  RETURN
END SUBROUTINE U_O_MVG


SUBROUTINE U_H_MVG(H, O, OR, OS, A, N)
!**************************************************************************************************
! 采用 Mualem (1976) 和 van Genuchten (1980) 公式计算土壤水势
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: H                         ! 土壤水势 (cm)

  ! 输入参数:
  REAL, INTENT(IN)  :: O                         ! 土壤含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  REAL, INTENT(IN)  :: N                         ! 土壤孔径分布指数

  ! 局部变量:
  REAL              :: M                         ! M = 1 - 1 / N
  REAL              :: SE                        ! 土壤有效饱和度

  H = 0.0
  IF (O < OS) THEN
    M = 1.0 - 1.0 / N
    SE = (O - OR) / (OS - OR)
    H = - 1.0 / A * (SE ** (-1.0 / M) - 1.0) ** (1.0 / N)
  END IF

  RETURN
END SUBROUTINE U_H_MVG


SUBROUTINE U_HYCAP_MVG(C, H, OR, OS, A, N)
!**************************************************************************************************
! 采用 Mualem (1976) 和 van Genuchten (1980) 公式计算土壤水容量
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: C                         ! 土壤水容量 (cm^-1)

  ! 输入参数:
  REAL, INTENT(IN)  :: H                         ! 土壤水势 (cm)
  REAL, INTENT(IN)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  REAL, INTENT(IN)  :: N                         ! 土壤孔径分布指数

  ! 局部变量:
  REAL              :: M                         ! M = 1 - 1 / N

  C = 0.0
  IF (H < 0.0) THEN
    M = 1.0 - 1.0 / N
    C = (OS - OR) * A * M * N * (- A * H) ** (N - 1.0) * &
        (1.0 + (- A * H) ** N) ** (- M - 1.0)
  END IF

  RETURN
END SUBROUTINE U_HYCAP_MVG


subroutine U_HYCON_BC(K, H, KS, A, N, L)
!**************************************************************************************************
! 采用 Brooks & Corey (1964) 模型计算土壤非饱和导水率
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: K                         ! 土壤非饱和导水率 (cm d^-1)

  ! 输入参数:
  real, intent(in)  :: H                         ! 土壤水势 (cm)
  real, intent(in)  :: KS                        ! 土壤饱和导水率 (cm d^-1)
  real, intent(in)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  real, intent(in)  :: N                         ! 土壤孔径分布指数
  real, intent(in)  :: L                         ! 土壤孔隙连结性参数

  ! 局部变量:
  real              :: SE                        ! 土壤有效饱和度

  K = KS
  if (H < (- 1.0 / A)) then
    SE = abs(A * H) ** (- N)
    K = KS * SE ** (2.0 / N + L + 2.0)
  end if

  return
end subroutine U_HYCON_BC


subroutine U_O_BC(O, H, OR, OS, A, N)
!**************************************************************************************************
! 采用 Brooks & Corey (1964) 模型计算土壤含水率
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: O                         ! 土壤含水率 (cm^3 cm^-3)

  ! 输入参数:
  real, intent(in)  :: H                         ! 土壤水势 (cm)
  real, intent(in)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  real, intent(in)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  real, intent(in)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  real, intent(in)  :: N                         ! 土壤孔径分布指数

  ! 局部变量:
  real              :: SE                        ! 土壤有效饱和度

  O = OS
  if (H < (- 1.0 / A)) then
    SE = abs(A * H) ** (- N)
    O = OR + SE * (OS - OR)
  end if

  return
end subroutine U_O_BC


subroutine U_H_BC(H, O, OR, OS, A, N)
!**************************************************************************************************
! 采用 Brooks & Corey (1964) 模型计算土壤水势
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: H                         ! 土壤水势 (cm)

  ! 输入参数:
  real, intent(in)  :: O                         ! 土壤含水率 (cm^3 cm^-3)
  real, intent(in)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  real, intent(in)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  real, intent(in)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  real, intent(in)  :: N                         ! 土壤孔径分布指数

  ! 局部变量:
  real              :: SE                        ! 土壤有效饱和度

  H = - 1.0 / A
  if (O < OS) then
    SE = (O - OR) / (OS - OR)
    H = - SE ** (- 1.0 / N) / A
  end if

  return
end subroutine U_H_BC


SUBROUTINE U_HYCAP_BC(C, H, OR, OS, A, N)
!**************************************************************************************************
! 采用 Brooks & Corey (1964) 模型计算土壤水容量
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: C                         ! 土壤水容量 (cm^-1)

  ! 输入参数:
  real, intent(in)  :: H                         ! 土壤水势 (cm)
  real, intent(in)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  real, intent(in)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  real, intent(in)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  real, intent(in)  :: N                         ! 土壤孔径分布指数

  !C = TINY(0.0)
  C = 0.0
  if (H < -1.0 / A) C = (OS - OR) * A * N * (- A * H) ** (- N - 1.0)

  return
end subroutine U_HYCAP_BC


subroutine U_Kos_HyCon(K, H, KS, A, N, L)
!**************************************************************************************************
! 采用 Kosugi (1996) 模型计算土壤非饱和导水率
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: K                         ! 土壤非饱和导水率 (cm d^-1)

  ! 输入参数:
  real, intent(in)  :: H                         ! 土壤水势 (cm)
  real, intent(in)  :: KS                        ! 土壤饱和导水率 (cm d^-1)
  real, intent(in)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  real, intent(in)  :: N                         ! 土壤孔径分布指数
  real, intent(in)  :: L                         ! 土壤孔隙连结性参数

  ! 局部变量:
  real              :: SE                        ! 土壤有效饱和度

  K = KS
  if (H < 0.0) then
    SE = erfc(log(h / A) / N)
    K = KS * SE ** L * (erfc(log(h / A) / N + N)) ** 2.0
  end if

  return
end subroutine U_Kos_HyCon


subroutine U_Kos_O(O, H, OR, OS, A, N)
!**************************************************************************************************
! 采用 Kosugi (1996) 模型计算计算土壤含水率
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: O                         ! 土壤含水率 (cm^3 cm^-3)

  ! 输入参数:
  real, intent(in)  :: H                         ! 土壤水势 (cm)
  real, intent(in)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  real, intent(in)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  real, intent(in)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  real, intent(in)  :: N                         ! 土壤孔径分布指数

  O = OS
  if (H < 0.0) then
    O = OR + (OS - OR) * erfc(log(h / A) / N)
  end if

  return
end subroutine U_Kos_O


subroutine U_Kos_H(H, O, OR, OS, A, N)
!**************************************************************************************************
! 采用 Kosugi (1996) 公式计算土壤水势
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: H                         ! 土壤水势 (cm)

  ! 输入参数:
  real, intent(in)  :: O                         ! 土壤含水率 (cm^3 cm^-3)
  real, intent(in)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  real, intent(in)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  real, intent(in)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  real, intent(in)  :: N                         ! 土壤孔径分布指数

  ! 局部变量:
  real, dimension(1) :: SE                        ! 土壤有效饱和度
  real               :: Y

  H = 0.0
  if (O < OS) then
    SE(1) = (O - OR) / (OS - OR)
    ! call vserfcinv(1, SE, Y)
    H = - A * exp(N * Y)
  end if

  return
end subroutine U_Kos_H


subroutine U_Kos_HyCap(C, H, OR, OS, A, N)
!**************************************************************************************************
! 采用 Kosugi (1996) 公式计算土壤水容量
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: C                         ! 土壤水容量 (cm^-1)

  ! 输入参数:
  real, intent(in)  :: H                         ! 土壤水势 (cm)
  real, intent(in)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  real, intent(in)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  real, intent(in)  :: A                         ! 土壤进气值的倒数 (cm^-1)
  real, intent(in)  :: N                         ! 土壤孔径分布指数

  ! 常数:
  real, parameter   :: PI = 4.0 * atan(1.0)      ! 圆周率

  ! 局部变量:
  real              :: M                         ! M = 1 - 1 / N

  C = tiny(0.0)
  if (H < 0.0) then
    C = (OR - OS) * 2.0 * exp(- log((abs(H) / A) ** 2.0) / N ** 2.0) / &
        abs(H) / N / sqrt(PI)
  end if

  return
end subroutine U_Kos_HyCap


SUBROUTINE U_RWUWSR_FED(WSR, H, T, HCRT, TCRT)
!**************************************************************************************************
! 采用 Feddes (1978) 公式计算根系吸水的水分胁迫响应函数值
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT)              :: WSR          ! 水分胁迫响应函数值

  ! 输入参数:
  REAL, INTENT(IN)               :: H            ! 土壤水势 (cm)
  REAL, INTENT(IN)               :: T            ! 潜在蒸腾速率 (cm d^-1)
  REAL, DIMENSION(7), INTENT(IN) :: HCRT         ! 根系吸水水分胁迫响应函数的临界水势 (cm)
  REAL, DIMENSION(7), INTENT(IN) :: TCRT         ! 根系吸水水分胁迫响应函数的临界蒸腾速率 (cm d^-1)

  ! 局部变量:
  REAL                           :: H1           ! 厌氧点的土壤水势 (cm)
  REAL                           :: H2           ! WSR = 1 的土壤水势上限
  REAL                           :: H3           ! WSR = 1 的土壤水势下限
  REAL                           :: H4           ! 凋萎点的土壤水势 (cm)
  REAL                           :: H3H          ! 与 T3H 对应的土壤水势
  REAL                           :: H3L          ! 与 T3L 对应的土壤水势
  REAL                           :: T3H          !
  REAL                           :: T3L          !

  H1  = HCRT(1)
  H2  = HCRT(2)
  H3H = HCRT(3)
  H3L = HCRT(4)
  H4  = HCRT(5)
  T3H = TCRT(6)
  T3L = TCRT(7)

  ! 计算 H3
  IF (T > T3H) THEN
    H3 = H3H
  ELSE IF (T > T3L) THEN
    H3 = H3H + (H3L - H3H) * (T3H - T) / (T3H - T3L)
  ELSE
    H3 = H3L
  END IF

  ! 计算根系吸水的水分胁迫系数
  IF (H < H4) THEN
    WSR = 0.0
  ELSE IF (H < H3) THEN
    WSR = (H - H4) / (H3 - H4)
  ELSE IF (H < H2) THEN
    WSR = 1.0
  ELSE IF (H < H1) THEN
    WSR = (H1 - H) / (H1 - H2)
  ELSE
    WSR = 0.0
  END IF

  RETURN
END SUBROUTINE U_RWUWSR_FED


SUBROUTINE U_KRWUS_VG1(WSR, H, HO, H50, P)
!**************************************************************************************************
! 采用 van Genuchten (1987) 公式计算根系吸水的水分胁迫响应函数值
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: WSR        ! 水分胁迫响应函数值

  ! 输入参数:
  REAL, INTENT(IN)  :: H          ! 土壤水势 (cm)
  REAL, INTENT(IN)  :: HO         ! Osmotic head (cm)
  REAL, INTENT(IN)  :: H50        ! Pressure head at which the water extraction rate is reduced 
                                  ! by 50% during conditions of negligible osmotic stress (cm)
  REAL, INTENT(IN)  :: P          ! Empirical parameter

  WSR = 1.0 / (1.0 + ((H + HO) / H50) ** P)

  RETURN
END SUBROUTINE U_KRWUS_VG1


SUBROUTINE U_KRWUS_VG2(WSR, H, HO, H50, HO50, P1, P2)
!**************************************************************************************************
! 采用 van Genuchten (1987) 公式计算根系吸水的水分胁迫响应函数值
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: WSR        ! 水分胁迫响应函数值

  ! 输入参数:
  REAL, INTENT(IN) :: H           ! 土壤水势 (cm)
  REAL, INTENT(IN) :: HO          ! Osmotic head (cm)
  REAL, INTENT(IN) :: H50         ! Pressure head at which the water extraction rate is reduced 
                                  ! by 50% during conditions of negligible osmotic stress (cm)
  REAL, INTENT(IN) :: HO50        ! Osmotic head at which the water extraction rate is reduced by
                                  ! 50% during conditions of negligible water stress (cm)
  REAL, INTENT(IN) :: P1          ! Empirical parameter
  REAL, INTENT(IN) :: P2          ! Empirical parameter

  WSR = 1.0 / (1.0 + (H / H50) ** P1) + 1.0 / (1.0 + (HO / HO50) ** P2)

  RETURN
END SUBROUTINE U_KRWUS_VG2


SUBROUTINE U_THCON_DV(THCON, O, OSAT, SAND, CLAY, THCONSAND, THCONCLAY, THCONOM, THCONWAT, THCONAIR)
!**************************************************************************************************
! 采用 De Vries (1963) 公式计算土壤导热率
!**************************************************************************************************
  implicit none

  ! 输出参数:
  REAL, INTENT(OUT) :: THCON                     ! 导热率 (J cm^-1 °C^-1 d^-1)

  ! 输入参数:
  REAL, INTENT(IN)  :: O                         ! 含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OSAT                      ! 饱和含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: SAND                      ! 砂粒占固相体积的比例
  REAL, INTENT(IN)  :: CLAY                      ! 粘粒占固相体积的比例
  REAL, INTENT(IN)  :: THCONSAND                 ! 砂粒的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, INTENT(IN)  :: THCONCLAY                 ! 粘粒的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, INTENT(IN)  :: THCONOM                   ! 有机质的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, INTENT(IN)  :: THCONWAT                  ! 水的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, INTENT(IN)  :: THCONAIR                  ! 空气的导热率 (J cm^-1 °C^-1 d^-1)

  ! 常数
  REAL, PARAMETER   :: XSANDWAT = 0.2461
  REAL, PARAMETER   :: XCLAYWAT = 0.7317
  REAL, PARAMETER   :: XOMWAT   = 1.2602
  REAL, PARAMETER   :: XWATWAT  = 1.0000
  REAL, PARAMETER   :: XAIRWAT  = 1.6062
  REAL, PARAMETER   :: XSANDAIR = 0.0143
  REAL, PARAMETER   :: XCLAYAIR = 0.6695
  REAL, PARAMETER   :: XOMAIR   = 0.4545
  REAL, PARAMETER   :: XWATAIR  = 0.1812
  REAL, PARAMETER   :: XAIRAIR  = 1.0000

  REAL, PARAMETER   :: ODRY = 0.02
  REAL, PARAMETER   :: OWET = 0.05

  ! 局部变量:
  REAL              :: FSAND                     ! 土壤砂粒体积含量
  REAL              :: FCLAY                     ! 土壤粘粒体积含量
  REAL              :: FOM                       ! 土壤有机质体积含量
  REAL              :: FAIR                      ! 土壤空气体积含量
  REAL              :: THCONWET
  REAL              :: THCONDRY

  FSAND = SAND * (1.0 - OSAT)
  FCLAY = CLAY * (1.0 - OSAT)
  FOM = (1.0 - OSAT) - FSAND - FCLAY
  FAIR = OSAT - O

  THCONWET = (XSANDWAT * FSAND * THCONSAND + XCLAYWAT * FCLAY * THCONCLAY + &
              XOMWAT * FOM * THCONOM + XWATWAT * O * THCONWAT + &
              XAIRWAT * FAIR * THCONAIR) / (XSANDWAT * FSAND + &
              XCLAYWAT * FCLAY + XOMWAT * FOM + XWATWAT * O + XAIRWAT * FAIR)

  THCONDRY = 1.25 * (XSANDAIR * FSAND * THCONSAND + XCLAYAIR * FCLAY * THCONCLAY + &
              XOMAIR * FOM * THCONOM + XWATAIR * O * THCONWAT + &
              XAIRAIR * FAIR * THCONAIR) / (XSANDAIR * FSAND + &
              XCLAYAIR * FCLAY + XOMAIR * FOM + XWATAIR * O + XAIRAIR * FAIR)

  THCON = THCONDRY + (THCONWET - THCONDRY) / (OWET - ODRY) * (O - ODRY)

  RETURN
END SUBROUTINE U_THCON_DV


SUBROUTINE U_THCON_CK(THCON, O, OR, OS, THCONSLD, THCONWAT, X, A, B)
!**************************************************************************************************
! 采用 Cote & Konrad (2005) 公式计算土壤导热率
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: THCON                     ! 土壤导热率 (J cm^-1 °C^-1 d^-1)

  ! 输入参数:
  REAL, INTENT(IN)  :: O                         ! 土壤含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OR                        ! 土壤残余含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OS                        ! 土壤饱和含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: THCONSLD                  ! 土壤固相的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, INTENT(IN)  :: THCONWAT                  ! 土壤水的导热率 (J cm^-1 °C^-1 d^-1)
  REAL, INTENT(IN)  :: X                         ! 与土壤类型有关的经验参数
  REAL, INTENT(IN)  :: A                         ! 土壤颗粒形状参数
  REAL, INTENT(IN)  :: B                         ! 土壤颗粒形状参数

  ! 局部变量:
  REAL              :: SE
  REAL              :: THCONR
  REAL              :: THCONDRY
  REAL              :: THCONSAT

  SE = (O - OR) / (OS - OR)
  THCONR = X * SE / (1.0 + (X - 1.0) * SE)
  THCONDRY = A * 10.0 ** (- B * OS)
  THCONSAT = THCONSLD ** (1.0 - OS) * THCONWAT ** OS
  THCON = THCONDRY + (THCONSAT - THCONDRY) * THCONR

  RETURN
END SUBROUTINE U_THCON_CK


SUBROUTINE U_THCON_LU(THCON, O, OSAT, SAND, CLAY, BD)
!**************************************************************************************************
! 采用 Lu (2014) 公式计算土壤导热率
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: THCON                     ! 土壤导热率 (J cm^-1 °C^-1 d^-1)

  ! 输入参数:
  REAL, INTENT(IN)  :: O                         ! 土壤含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OSAT                      ! 土壤饱和含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: SAND                      ! 土壤固相的砂粒体积含量
  REAL, INTENT(IN)  :: CLAY                      ! 土壤固相的粘粒体积含量
  REAL, INTENT(IN)  :: BD                        ! 土壤干容重 (g cm^-3)

  ! 局部变量:
  REAL              :: BET                       ! 经验参数
  REAL              :: ALP                       ! 经验参数

  BET = 1.97 * SAND + 1.87 * BD - 1.36 * SAND * BD - 0.95
  ALP = 0.67 * CLAY + 0.24
  THCON = - 0.56 * OSAT + 0.51 + EXP(BET - O ** (- ALP))

  ! W 换算为 J d^-1
  THCON = THCON * 864.0 * 274.15

  RETURN
END SUBROUTINE U_THCON_LU


SUBROUTINE U_THCAP_DV(THCAP, O, OSAT, SAND, CLAY)
!**************************************************************************************************
! 采用 De Vries (1963) 公式计算热容量
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL, INTENT(OUT) :: THCAP                     ! 土壤热容量 (J cm^-3 °C^-1)

  ! 输入参数:
  REAL, INTENT(IN)  :: O                         ! 土壤含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: OSAT                      ! 土壤饱和含水率 (cm^3 cm^-3)
  REAL, INTENT(IN)  :: SAND                      ! 土壤固相的砂粒体积含量
  REAL, INTENT(IN)  :: CLAY                      ! 土壤固相的粘粒体积含量

  ! 常数
  REAL, PARAMETER   :: THCAPSAND = 2.128         ! 土壤砂粒的热容量 (J cm^-1 °C^-1 d^-1)
  REAL, PARAMETER   :: THCAPCLAY = 2.385         ! 土壤粘粒的热容量 (J cm^-1 °C^-1 d^-1)
  REAL, PARAMETER   :: THCAPOM   = 2.496         ! 土壤有机质的热容量 (J cm^-1 °C^-1 d^-1)
  REAL, PARAMETER   :: THCAPWAT  = 4.180         ! 土壤水的热容量 (J cm^-1 °C^-1 d^-1)
  REAL, PARAMETER   :: THCAPAIR  = 0.001212      ! 土壤空气的热容量 (J cm^-1 °C^-1 d^-1)

  ! 局部变量:
  REAL              :: FSAND                     ! 土壤砂粒体积含量
  REAL              :: FCLAY                     ! 土壤粘粒体积含量
  REAL              :: FOM                       ! 土壤有机质体积含量
  REAL              :: FAIR                      ! 土壤空气体积含量

  FSAND = SAND * (1.0 - OSAT)
  FCLAY = CLAY * (1.0 - OSAT)
  FOM = (1.0 - OSAT) - FSAND - FCLAY
  FAIR = OSAT - O

  THCAP = FSAND * THCAPSAND + FCLAY * THCAPCLAY + FOM * THCAPOM + &
          O * THCAPWAT + FAIR * THCAPAIR

  RETURN
END SUBROUTINE U_THCAP_DV


FUNCTION SVP(AT)
!**************************************************************************************************
! 计算给定气温的饱和水汽压
!**************************************************************************************************
  IMPLICIT NONE

  ! 输出参数:
  REAL             :: SVP                        ! 饱和水汽压 (kPa)

  ! 输入参数:
  REAL, INTENT(IN) :: AT                         ! 气温 (°C)

  SVP = 0.6108 * EXP(17.27 * AT / (AT + 237.3))

  RETURN
END FUNCTION SVP


subroutine U_RadExt(RA, DHMax, Lat, DOY)
!**************************************************************************************************
! 计算天文辐射量和最大可能日照时数
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real,    intent(out) :: RA                     ! 天文辐射量 (MJ m^-2 d^-1)
  real,    intent(out) :: DHMax                  ! 最大可能日照时数 (h)

  ! 输入参数:
  real,    intent(IN)  :: Lat                    ! 纬度 (°)
  integer, intent(IN)  :: DOY                    ! 年日序 (d)

  ! 常数:
  real, parameter      :: PI = 4.0 * atan(1.0)   ! 圆周率

  ! 局部变量:
  real                 :: Y                      ! 纬度, 以弧度表示 (rad)
  real                 :: J                      !
  real                 :: DR                     ! Inverse relative distance Earth-Sun
  real                 :: D                      ! Solar decimation (rad)
  real                 :: W                      ! Sunset hour angle (rad)

  Y = LAT * PI / 180.0
  J = 2.0 * PI * DOY / 365.0
  DR = 1.0 + 0.033 * cos(J)
  D = 0.409 * sin(J - 1.39)
  W = acos(-tan(Y) * tan(D))
  RA = 118.08 / PI * DR * (W * sin(Y) * sin(D) + sin(W) * cos(Y) * cos(D))
  DHMax = 24.0 / PI * W

  return
end subroutine U_RadExt


subroutine U_Rad(Rad, DH, DHMax, RA, KAngstrom)
!**************************************************************************************************
! 计算太阳辐射量和最大可能日照时数
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: Rad                       ! 太阳辐射量 (MJ m^-2 d^-1)

  ! 输入参数:
  real, intent(in)  :: DH                        ! 日照时数 (h)
  real, intent(in)  :: DHMax                     ! 最大可能日照时数 (h)
  real, intent(in)  :: RA                        ! 天文辐射 (MJ m^-2 d^-1)
  real, dimension(2), intent(in)  :: KAngstrom   ! Angstrom-Prescott 公式系数

  Rad = (KAngstrom(1) + KAngstrom(2) * DH / DHMax) * RA

  return
end subroutine U_Rad


subroutine U_RadNet(RN, ATMax, ATMin, RHMean, Rad, RA, KAngstrom, Albedo)
!**************************************************************************************************
! 计算净辐射量
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real,               intent(out) :: RN          ! 净辐射量 (MJ m^-2 d^-1)

  ! 输入参数:
  real,               intent(in)  :: ATMax       ! 最高气温 (°C)
  real,               intent(in)  :: ATMin       ! 最低气温 (°C)
  real,               intent(in)  :: RHMean      ! 平均相对湿度 (%)
  real,               intent(in)  :: Rad         ! 太阳辐射量 (MJ m^-2 d^-1)
  real,               intent(in)  :: RA          ! 天文辐射量 (MJ m^-2 d^-1)
  real, dimension(2), intent(in)  :: KAngstrom   ! Angstrom-Prescott 公式系数
  real,               intent(in)  :: Albedo      ! 地表反射率

  ! 函数声明:
  real              :: SVP

  ! 常数:
  real, parameter   :: SB = 0.000000004903       ! Stefan-Boltzmann 常数

  ! 局部变量:
  real              :: RSO                       ! 晴空辐射量 (MJ m^-2 d^-1)
  real              :: RNS                       ! 净短波辐射量 (MJ m^-2 d^-1)
  real              :: RNL                       ! 净长波辐射量 (MJ m^-2 d^-1)
  real              :: ES                        ! 饱和水汽压 (kPa)
  real              :: EA                        ! 实际水汽压 (kPa)

  ES = (SVP(ATMax) + SVP(ATMin)) * 0.5
  EA = ES * RHMean * 0.01
  RSO = (KAngstrom(1) + KAngstrom(2)) * RA
  RNS = (1.0 - Albedo) * Rad
  RNL = SB * ((ATMax + 273.16) ** 4.0 + (ATMin + 273.16) ** 4.0) * 0.5 * &
       (0.34 - 0.14 * sqrt(EA)) * (1.35 * Rad / RSO - 0.35)
  RN = RNS - RNL

  return
end subroutine


subroutine U_WindSpeed2(WS2, WS, H)
!**************************************************************************************************
! 计算地面 2m 高处风速
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: WS2                       ! 地面 2 m 高处风速 (m s^-1)

  ! 输入参数:
  real, intent(in)  :: H                         ! 风速测量高度 (m)
  real, intent(in)  :: WS                        ! 风速 (m s^-1)

  WS2 = WS
  if (H > 2.0) WS2 = WS * 4.87 / log(67.8 * H - 5.42)

  return
end subroutine


subroutine U_FaoPenMon(ETRef, Alt, ATMax, ATMin, RHMean, WS2, RN)
!**************************************************************************************************
! 采用 FAO Penman-Monteith 公式计算参考作物蒸发蒸腾量
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: ETRef                     ! 参考作物蒸发蒸腾量 (mm d^-1)

  ! 输入参数:
  real, intent(in)  :: Alt                       ! 海拔 (m)
  real, intent(in)  :: ATMax                     ! 最高气温 (°C)
  real, intent(in)  :: ATMin                     ! 最低气温 (°C)
  real, intent(in)  :: RHMean                    ! 平均相对湿度 (%)
  real, intent(in)  :: WS2                       ! 地面 2m 高处风速 (m s^-1)
  real, intent(in)  :: RN                        ! 日净辐射量 (MJ m^-2 d^-1)

  ! 函数声明:
  real              :: SVP

  ! 局部变量:
  real              :: ATMean                    ! 平均气温 (°C)
  real              :: Gam                       ! Psychrometric constant (kPa °C^-1)
  real              :: Del                       ! Slope vapour pressure curve (kPa °C^-1)
  real              :: ES                        ! 饱和水汽压 (kPa)
  real              :: EA                        ! 实际水汽压 (kPa)
  real              :: ETRad                     ! Radiation term (mm d^-1)
  real              :: ETAero                    ! Aerodynamic term (mm d^-1)

  ATMean = (ATMin + ATMax) * 0.5
  Gam = 0.000665 * 101.3 * (1.0 - 0.0065 * Alt / 293.0) ** 5.26
  Del = 4098.0 * SVP(ATMean) / (ATMean + 237.3) ** 2.0
  ES = (SVP(ATMin) + SVP(ATMax)) * 0.5
  EA = ES * RHMean * 0.01
  ETRad = 0.408 * Del * RN / (Del + Gam * (1.0 + 0.34 * WS2))
  ETAero = Gam * 900.0 / (ATMean + 273.0) * WS2 * (ES - EA) / (Del + Gam * (1.0 + 0.34 * WS2))
  ETRef = ETRad + ETAero
  if(ETRef < 0.0) ETRef = 0.0

  ! 单位由 mm 转化为 cm
  ETRef = ETRef * 0.1

  return
end subroutine U_FaoPenMon


subroutine U_FaoHarg(ETRef, ATMax, ATMin, RA)
!**************************************************************************************************
! 采用 FAO Hargreaves 公式计算参考作物蒸发蒸腾量
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: ETRef                     ! 参考作物蒸发蒸腾量 (mm d^-1)

  ! 输入参数:
  real, intent(in)  :: ATMax                     ! 最高气温 (°C)
  real, intent(in)  :: ATMin                     ! 最低气温 (°C)
  real, intent(in)  :: RA                        ! 日天文辐射量 (MJ m^-2 d^-1)

  ! 局部变量:
  real              :: ATMean                    ! 平均气温 (°C)
  real              :: Lam                       ! 气化潜热 (MJ kg^-1)

  ATMean = (ATMin + ATMax) * 0.5
  Lam = 2.501 - 0.002361 * ATMean
  ETRef = 0.0023 * (ATMean + 17.8) * sqrt(ATMax - ATMin) * RA

  ! 单位由 mm 转化为 cm
  ETRef = ETRef * 0.1

  return
end subroutine U_FaoHarg


subroutine U_PenMon(ETPot, Alt, ATMax, ATMin, RHMean, WS2, RN, G, Res1, Res2)
!**************************************************************************************************
! 采用 Penman-Monteth (1981) 公式计算潜在蒸发蒸腾量
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: ETPot                     ! 潜在蒸发蒸腾量 (mm d^-1)

  ! 输入参数:
  real, intent(in)  :: Alt                       ! 海拔 (m)
  real, intent(in)  :: ATMax                     ! 最高气温 (°C)
  real, intent(in)  :: ATMin                     ! 最低气温 (°C)
  real, intent(in)  :: RHMean                    ! 平均相对湿度 (%)
  real, intent(in)  :: RN                        ! 净辐射量 (MJ m^-2 d^-1)
  real, intent(in)  :: WS2                       ! 地面 2m 高处风速 (m s^-1)
  real, intent(in)  :: G                         ! 土壤热通量 (MJ m^-2 d^-1)
  real, optional    :: Res1                      ! 冠层阻力 (s m^-1)
  real, optional    :: Res2                      ! 空气动力学阻力 (s m^-1)

  ! 函数声明:
  real              :: SVP

  ! 局部变量:
  real              :: ATMean                    ! 平均气温 (°C)
  real              :: P                         ! 大气压 (kPa)
  real              :: Gam                       ! 湿度计常数 (kPa °C^-1)
  real              :: Del                       ! 饱和水汽压曲线的斜率 (kPa °C^-1)
  real              :: Lam                       ! 气化潜热 (MJ kg^-1)
  real              :: Rho                       ! 空气密度 (kg m^-3)
  real              :: ES                        ! 饱和水汽压 (kPa)
  real              :: EA                        ! 实际水汽压 (kPa)

  if (.not. present(Res1)) Res1 = 70.0
  if (.not. present(Res2)) Res2 = 208.0 / WS2

  ATMean = (ATMin + ATMax) * 0.5
  P = 101.3 * (1.0 - 0.0065 * Alt / 293.0) ** 5.26
  Gam = 0.000665 * P
  Rho = 3.486 * P / (1.01 * (ATMean + 273))
  Del = 4098.0 * SVP(ATMean) / (ATMean + 237.3) ** 2.0
  Lam = 2.501 - 0.002361 * ATMean
  ES = (SVP(ATMin) + SVP(ATMax)) * 0.5
  EA = ES * RHMean * 0.01
  ETPot = (Del * (RN - G) + 0.001013 * Rho * (ES - EA) / Res2) / &
          (Del + Gam * (1.0 + Res1 / Res2)) / Lam

  ! 单位由 mm 转化为 cm
  ETPot = ETPot * 0.1

  return
end subroutine U_PenMon


subroutine U_PriTay(ETPot, Alt, ATMax, ATMin, RHMean, RN, G, Alp)
!**************************************************************************************************
! 采用 Priestley-Taylor (1972) 公式计算潜在蒸发蒸腾量
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: ETPot                     ! 潜在蒸发蒸腾量 (mm d^-1)

  ! 输入参数:
  real, intent(in)  :: Alt                       ! 海拔 (m)
  real, intent(in)  :: ATMax                     ! 最高气温 (°C)
  real, intent(in)  :: ATMin                     ! 最低气温 (°C)
  real, intent(in)  :: RHMean                    ! 平均相对湿度 (%)
  real, intent(in)  :: RN                        ! 净辐射量 (MJ m^-2 d^-1)
  real, intent(in)  :: G                         ! 土壤热通量 (MJ m^-2 d^-1)
  real, optional    :: Alp                       ! 调整系数

  ! 函数声明:
  real              :: SVP

  ! 局部变量:
  real              :: ATMean                    ! 平均气温 (°C)
  real              :: P                         ! 大气压 (kPa)
  real              :: Gam                       ! 湿度计常数 (kPa °C^-1)
  real              :: Del                       ! 饱和水汽压曲线的斜率 (kPa °C^-1)
  real              :: Rho                       ! 空气密度 (kg m^-3)
  real              :: Lam                       ! 气化潜热 (MJ kg^-1)

  if (.not. present(Alp)) Alp = 1.25
  ATMean = (ATMin + ATMax) * 0.5
  P = 101.3 * (1.0 - 0.0065 * Alt / 293.0) ** 5.26
  Gam = 0.000665 * P
  Del = 4098.0 * SVP(ATMean) / (ATMean + 237.3) ** 2.0
  Rho = 3.486 * P / (1.01 * (ATMean + 273))
  Lam = 2.501 - 0.002361 * ATMean
  ETPot = Alp * Del / (Del + Gam) * (RN - G) / Lam

  ! 单位由 mm 转化为 cm
  ETPot = ETPot * 0.1

  return
end subroutine U_PriTay


subroutine U_ResAero(Res, H, WS, WH, HH)
!**************************************************************************************************
! 计算空气动力学阻力
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: Res                       ! 空气动力学阻力 (s m^-1)

  ! 输入参数:
  real, intent(in)  :: H                         ! 作物高度 (cm)
  real, intent(in)  :: WS                        ! 风速 (m s^-1)
  real, intent(in)  :: WH                        ! 风速测量高度 (cm)
  real, intent(in)  :: HH                        ! 湿度测量高度 (cm)

  ! 局部变量
  real              :: D                         ! Zero plane displacement height (cm)
  real              :: Zom                       ! Roughness length governing momentum transfer (cm)
  real              :: Zoh                       ! Roughness length governing transfer of heat and vapour (cm)

  D = 2.0 / 3.0 * H
  Zom = 0.123 * H
  Zoh = 0.1 * Zom
  Res = log((WH - D) / Zom)* log((HH - D) / Zoh) / 0.1681 / WS

 return
end subroutine U_ResAero


subroutine U_ResSurf(Res, ResL, LAI)
!**************************************************************************************************
! 计算表面阻力
!**************************************************************************************************
  implicit none

  ! 输出参数:
  real, intent(out) :: Res                       ! 表面阻力 (mm d^-1)

  ! 输入参数:
  real, intent(in)  :: ResL                      ! Bulk stomatal resistance of the well-illuminated leaf (s m-1)
  real, intent(in)  :: LAI                       ! 叶面积指数

  Res = ResL / LAI / 2.0

 return
end subroutine U_ResSurf

