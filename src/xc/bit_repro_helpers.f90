!> The arithmetic machinery behind the bit_repro public functions: exact-rounding
!! primitives, the double-double log2/exp2 core, binary exponentiation, and the
!! cuberoot rescalers.  Bodies are verbatim from the pre-submodule layout; the
!! split is bit-neutral and gated by the harness checksums.
submodule (bit_repro) bit_repro_helpers

implicit none

! IEEE-754 binary64 layout constants
integer, parameter :: fraclen = 52  !< Bits in the fraction (mantissa)
integer, parameter :: explen  = 11  !< Bits in the exponent
integer, parameter :: signbit = 63  !< Position of the sign bit
integer, parameter :: expbit  = fraclen !< Position of the lowest exponent bit
integer(kind=int64), parameter :: bias = 1023_int64 !< Exponent bias

contains

!> exponent()/fraction() split that is well-defined for denormal x > 0 on every
!! compiler: denormals are pre-scaled by 2**54 (an exact multiply) first, because
!! the intrinsics' behavior on denormals is processor-dependent (and wrong under
!! nvfortran).  For normal x this is bit-identical to the plain intrinsics.
elemental module subroutine norm_split(x, m, k)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in)  :: x  !< The argument, x > 0, finite
  real(wp), intent(out) :: m  !< fraction(x), in [0.5, 1)
  integer, intent(out) :: k !< exponent(x), corrected for denormals
  real(wp), parameter :: two_p54 = 18014398509481984.0_wp  ! 2**54, exact
  real(wp) :: xn
  integer :: kb
  xn = x ; kb = 0
  if (xn < tiny(xn)) then ; xn = xn * two_p54 ; kb = -54 ; endif
  k = exponent(xn) + kb ; m = fraction(xn)
end subroutine norm_split

!> scale(v, n) with reproducible behavior at BOTH edges of the exponent range:
!! Fortran leaves SCALE processor-dependent outside the model, and nvfortran
!! (a) flushes denormal results to zero where ifx rounds them, and (b) returns
!! Inf for n > 1023 even when v < 1 makes the result representable (measured:
!! exp(709.774) -> Inf under nv, correct 1.78e308 under ifx).  The last 64
!! doublings/halvings are therefore done with one IEEE multiply by an exact
!! power of two, which every conforming platform rounds -- and overflows --
!! identically.
elemental function scale_reprod(v, n) result(sv)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in) :: v     !< The value to be scaled, |v| in ~[0.5, 2.5)
  integer, intent(in) :: n  !< The power of two to scale by
  real(wp) :: sv                !< v * 2**n, correctly rounded incl. denormals/overflow
  real(wp), parameter :: two_m64 = 2.0_wp**(-64)  ! exact
  real(wp), parameter :: two_p64 = 2.0_wp**64     ! exact
  if (n > 960) then
    sv = scale(v, n - 64) * two_p64
  elseif (n >= -960) then
    sv = scale(v, n)
  else
    sv = scale(v, n + 64) * two_m64
  endif
end function scale_reprod

!> x**n for n >= 0 by binary exponentiation: pure IEEE multiplies, so it is
!! reproducible everywhere, and ipow(x,2) == RN(x*x) exactly.
elemental module function ipow(x, n) result(p)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in) :: x    !< The base
  integer, intent(in) :: n !< The non-negative exponent
  real(wp) :: p                !< x**|n|, or its reciprocal for n < 0

  real(wp) :: b    ! Running square
  integer :: m ! Remaining exponent bits

  m = abs(n)
  if (m == 0) then
    p = 1.0_wp
  else
    p = 1.0_wp ; b = x
    do while (m > 1)
      if (btest(m, 0)) p = p * b
      b = b * b
      m = shiftr(m, 1)
    enddo
    p = p * b
  endif
  if (n < 0) p = 1.0_wp / p
end function ipow


!> Knuth two_sum: s + e == a + b exactly, branch-free.
elemental module subroutine two_sum(a, b, s, e)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in)  :: a, b
  real(wp), intent(out) :: s, e
  real(wp) :: bb
  s = a + b
  bb = s - a
  e = (a - (s - bb)) + (b - bb)
end subroutine two_sum

!> Dekker/Veltkamp two_prod WITHOUT fma: p + e == a*b exactly.
!! Requires the compiler not to contract these expressions into fma (build with
!! -fno-fma / -Mnofma / -gpu=nofma) -- contraction would still be exact here but
!! would change the split arithmetic bit-for-bit between compilers.
elemental module subroutine two_prod(a, b, p, e)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in)  :: a, b
  real(wp), intent(out) :: p, e
  real(wp), parameter :: splitter = 134217729.0_wp  ! 2**27 + 1
  real(wp) :: a1, a2, b1, b2, ta, tb
  p = a * b
  ta = splitter * a
  a1 = ta - (ta - a) ; a2 = a - a1
  tb = splitter * b
  b1 = tb - (tb - b) ; b2 = b - b1
  e = ((a1*b1 - p) + a1*b2 + a2*b1) + a2*b2
end subroutine two_prod

!> log2(x) for finite normal/denormal x > 0 as a double-double pair (h, l),
!! accurate to ~2**-85 relative: exponent/fraction split, mantissa recentered to
!! [sqrt(1/2), sqrt(2)), atanh series with the leading term carried in dd.
elemental module subroutine log2_dd(x, h, l)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in)  :: x  !< The argument, x > 0, finite
  real(wp), intent(out) :: h  !< High part of log2(x)
  real(wp), intent(out) :: l  !< Low part of log2(x)

  real(wp), parameter :: sqrt2_2 = 0.70710678118654752440_wp ! sqrt(1/2)
  ! 2/ln(2) as a double-double constant
  real(wp), parameter :: two_invln2_hi = 2.8853900817779268_wp
  real(wp), parameter :: two_invln2_lo = 4.0710547481862066e-17_wp
  ! Reciprocal odd integers for the atanh series
  real(wp), parameter :: a3=1.0_wp/3.0_wp,   a5=1.0_wp/5.0_wp,   a7=1.0_wp/7.0_wp,   a9=1.0_wp/9.0_wp
  real(wp), parameter :: a11=1.0_wp/11.0_wp, a13=1.0_wp/13.0_wp, a15=1.0_wp/15.0_wp, a17=1.0_wp/17.0_wp
  real(wp), parameter :: a19=1.0_wp/19.0_wp, a21=1.0_wp/21.0_wp
  real(wp) :: m          ! Mantissa recentered to [sqrt(1/2), sqrt(2))
  real(wp) :: u          ! m - 1, exact by Sterbenz
  real(wp) :: vh, vl     ! m + 1 as a dd pair (exact via two_sum)
  real(wp) :: s, sl      ! (m-1)/(m+1) as a dd pair
  real(wp) :: ph, pe     ! two_prod scratch
  real(wp) :: lh2, ll2   ! Leading term (2/ln2)*s as a dd pair
  real(wp) :: rv         ! 1/vh, so the routine needs only one division
  real(wp) :: s2, s4, s8 ! Powers of s for the Estrin evaluation
  real(wp) :: tail       ! The series tail (double precision suffices)
  real(wp) :: h1, l1     ! Intermediate dd sum
  integer :: k       ! Binary exponent of x

  call norm_split(x, m, k)
  if (m < sqrt2_2) then ; m = m + m ; k = k - 1 ; endif

  u = m - 1.0_wp                       ! exact: m in [0.70, 1.42)
  call two_sum(m, 1.0_wp, vh, vl)      ! v = m + 1 exactly as (vh, vl)
  rv = 1.0_wp / vh                     ! one division for the whole routine
  s = u * rv                        ! quotient to ~1 ulp ...
  call two_prod(s, vh, ph, pe)
  sl = ((u - ph) - pe - s*vl) * rv  ! ... refined: (s, sl) = u/v to ~2**-105

  ! Leading term L = (2/ln2) * s in dd
  call two_prod(s, two_invln2_hi, lh2, ll2)
  ll2 = ll2 + (s*two_invln2_lo + sl*two_invln2_hi)

  ! Series tail: L * s2*P(s2), |tail| <= 0.005, double is ample.  P is evaluated
  ! in a FIXED Estrin grouping (explicit parentheses; this order is part of the
  ! reproducibility contract) to shorten the dependency chain.
  s2 = s*s ; s4 = s2*s2 ; s8 = s4*s4
  tail = lh2 * (s2 * ( ((a3 + s2*a5) + s4*(a7 + s2*a9)) + &
                       s8*(((a11 + s2*a13) + s4*(a15 + s2*a17)) + s8*(a19 + s2*a21)) ))

  call two_sum(lh2, tail, h1, l1)
  l1 = l1 + ll2
  call two_sum(real(k, wp), h1, h, l)   ! add the exponent, keeping the residual
  l = l + l1
end subroutine log2_dd

!> 2**(h+l) for a double-double exponent, |l| << 1.  Table-based reduction:
!! h = n + i/32 + r with |r| <= 1/64, so 2**h = 2**n * T(i) * 2**r with T(i) a
!! 32-entry double-double table and 2**r a short degree-6 polynomial.  The
!! T(i)*poly product is carried through two_prod, keeping the assembly error
!! near half an ulp.  Overflow -> Inf, underflow -> 0/denormal via scale().
elemental module function exp2_pair(h, l) result(e2)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in) :: h  !< High part of the exponent
  real(wp), intent(in) :: l  !< Low part of the exponent
  real(wp) :: e2             !< 2**(h+l)

  ! T(i) = 2**(i/32) as double-double (generated at 60-digit precision)
  real(wp), parameter :: t2hi(0:31) = [ &
    1.0_wp, 1.0218971486541166_wp, 1.0442737824274138_wp, 1.0671404006768237_wp, &
    1.0905077326652577_wp, 1.1143867425958924_wp, 1.1387886347566916_wp, 1.1637248587775775_wp, &
    1.189207115002721_wp, 1.215247359980469_wp, 1.241857812073484_wp, 1.2690509571917332_wp, &
    1.2968395546510096_wp, 1.3252366431597413_wp, 1.3542555469368927_wp, 1.383909881963832_wp, &
    1.4142135623730951_wp, 1.4451808069770467_wp, 1.4768261459394993_wp, 1.5091644275934228_wp, &
    1.5422108254079407_wp, 1.5759808451078865_wp, 1.6104903319492543_wp, 1.645755478153965_wp, &
    1.681792830507429_wp, 1.718619298122478_wp, 1.7562521603732995_wp, 1.7947090750031072_wp, &
    1.8340080864093424_wp, 1.8741676341103_wp, 1.9152065613971474_wp, 1.9571441241754002_wp ]
  real(wp), parameter :: t2lo(0:31) = [ &
    0.0_wp, 5.109225028973444e-17_wp, 8.551889705537965e-17_wp, -7.899853966841582e-17_wp, &
    -3.046782079812471e-17_wp, 1.0410278456845571e-16_wp, 8.912812676025408e-17_wp, 3.8292048369240935e-17_wp, &
    3.982015231465646e-17_wp, -7.712630692681488e-17_wp, 4.658027591836937e-17_wp, 2.667932131342186e-18_wp, &
    2.5382502794888315e-17_wp, -2.8587312100388614e-17_wp, 7.70094837980299e-17_wp, -6.770511658794786e-17_wp, &
    -9.667293313452913e-17_wp, -3.0237581349939873e-17_wp, -3.483994556892796e-17_wp, -1.016455327754295e-16_wp, &
    7.949834809697621e-17_wp, -1.0136916471278304e-17_wp, 2.4707192569797888e-17_wp, -1.0125679913674773e-16_wp, &
    8.199010020581497e-17_wp, -1.851380418263111e-17_wp, 2.960140695448873e-17_wp, 1.8227458427912087e-17_wp, &
    3.283107224245627e-17_wp, -6.122763413004143e-17_wp, -1.0619946056195963e-16_wp, 8.960767791036668e-17_wp ]
  ! Coefficients ln2**j / j! for 2**r over |r| <= 1/64
  real(wp), parameter :: e2c1 = 0.6931471805599453_wp,    e2c2 = 0.24022650695910072_wp
  real(wp), parameter :: e2c3 = 0.05550410866482158_wp,   e2c4 = 0.009618129107628477_wp
  real(wp), parameter :: e2c5 = 0.0013333558146428443_wp, e2c6 = 0.0001540353039338161_wp
  real(wp) :: w    ! h*32
  real(wp) :: r    ! Fractional remainder of the exponent, |r| <= 1/64
  real(wp) :: q    ! 2**r - 1 from the polynomial, |q| <= 0.011
  real(wp) :: v    ! The assembled product T * 2**r before the l correction
  integer :: n32, idx, n

  if (h > 1100.0_wp) then
    e2 = huge(h)                     ! +Inf via runtime IEEE overflow (a constant
    e2 = e2 * 2.0_wp                    ! expression here is a compile error on gfortran)
  elseif (h < -1130.0_wp) then
    e2 = 0.0_wp
  else
    w = h * 32.0_wp
    n32 = nint(w)
    r = h - real(n32, wp) * 0.03125_wp      ! exact: n32/32 is an exact multiple of 2**-5
    idx = iand(n32, 31)              ! n32 = 32*n + idx with 0 <= idx < 32 for any
    n = (n32 - idx) / 32             ! sign; the division is exact (nvfortran has no shifta)
    ! q = 2**r - 1, never materializing the 1, assembled as
    !   v = T + (T*q + Tlo);  result = v + v*(l*ln2)
    ! The l correction must scale the FULL product v = T*2**r: applying it to T
    ! alone drops the q*l*ln2 cross term, which reaches ~2 ulp where |l| is
    ! largest (|h| in [256,1024), measured as the 3-ulp exp tail).
    q = r*(e2c1 + r*(e2c2 + r*(e2c3 + r*(e2c4 + r*(e2c5 + r*e2c6)))))
    v = t2hi(idx) + (t2hi(idx)*q + t2lo(idx))
    e2 = scale_reprod(v + v*(l*e2c1), n)
  endif
end function exp2_pair

!> Rescale `a` to the range [0.125, 1) while computing its cube-root exponent and
!! sign, by direct manipulation of the IEEE-754 bit representation.
pure module subroutine rescale_cbrt(a, x, e_r, s_a)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in) :: a  !< The number to be rescaled for cube-root computation
  real(wp), intent(out) :: x !< The rescaled value of `a` in the range [0.125, 1)
  integer(kind=int64), intent(out) :: e_r !< The integral component of the cube-root exponent of `a`
  integer(kind=int64), intent(out) :: s_a !< Sign bit of `a`; nonzero indicates negative

  integer(kind=int64) :: xb  ! Bit representation of `a`
  integer(kind=int64) :: e_a ! Exponent of `a`
  integer(kind=int64) :: e_x ! Exponent of `x`

  xb = transfer(a, 1_int64)
  s_a = ibits(xb, signbit, 1)
  e_a = ibits(xb, expbit, explen) - bias
  ! e = 3*(floor(e/3)+1) + (modulo(e,3) - 3); the last term is the exponent of x.
  e_r = (e_a + sign(1_int64, e_a) + 2) / 3
  e_x = e_a - e_r * 3
  call mvbits(e_x + bias, 0, explen + 1, xb, fraclen)
  x = transfer(xb, 1._wp)
end subroutine rescale_cbrt

!> Undo the rescaling of a real number back to its original base.
pure module function descale(x, e_a, s_a) result(a)
  !$omp declare target
  !$acc routine seq
  real(wp), intent(in) :: x !< The rescaled value which is to be restored
  integer(kind=int64), intent(in) :: e_a !< Exponent of the unscaled value
  integer(kind=int64), intent(in) :: s_a !< Sign bit of the unscaled value
  real(wp) :: a !< Restored value with the corrected exponent and sign

  integer(kind=int64) :: xb  ! Bit representation
  integer(kind=int64) :: e_x ! Biased exponent of x

  xb = transfer(x, 1_int64)
  e_x = ibits(xb, expbit, explen)
  call mvbits(e_a + e_x, 0, explen, xb, expbit)
  call mvbits(s_a, 0, 1, xb, signbit)
  a = transfer(xb, 1._wp)
end function descale

end submodule bit_repro_helpers
