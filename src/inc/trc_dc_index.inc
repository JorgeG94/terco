! TRC_DC_INDEX(i): how a `do concurrent` body hands its loop index (or a
! scalar it computed) to an item routine.
!
! ifx 2025.0 and 2026.0 at -O2 and above treat a call that takes the loop
! index by reference as loop-invariant when they do not inline the callee:
! every element gets the value of the first iteration (RI-CCSD T1 wrong by
! 3e-2 relative, (T) by 1e-2 Eh in the host build; -O1 and ifx 2022.1 are
! right). Passing an expression, (i), cures it. Every other compiler gets
! the plain variable, as written everywhere else in terco: nvfortran is the
! compiler this code is built for, and its path does not change for a defect
! of another compiler. (A VALUE dummy also cures ifx, but nvfortran 26.5 then
! gives wrong amplitudes in the host build.) Remove once ifx is fixed: the
! reproducer is in the "submit to Intel" issue.
#ifdef __INTEL_LLVM_COMPILER
#define TRC_DC_INDEX(i) (i)
#else
#define TRC_DC_INDEX(i) i
#endif
