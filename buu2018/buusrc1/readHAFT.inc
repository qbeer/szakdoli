c  HAFT declaration of acceptance matrix arrays
c
c  Code version 1.0 of February 10, 2006
c
c
c  The dimensions MUST match all array sizes in the file!
c
      integer size, nids
      parameter (size=250000) ! <<== change if < xdim*ydim*zdim
      parameter (nids=14)     ! <<== change if < id
c
      character*80 fname, fname2
      character*80 comment, comment2
      integer xdim(nids), ydim(nids), zdim(nids)
      real*8 dp(nids), dth(nids), dph(nids)
      real*8 pmin(nids), pmax(nids), thmin(nids),
     +       thmax(nids), phmin(nids), phmax(nids)
      real*8 matrix2(size) ! matrices are declared for
      real*8 matrix3(size) ! e+, e-, pi+, pi-, K+, K- and p
      real*8 matrix8(size)
      real*8 matrix9(size)
      real*8 matrix10(size)
      real*8 matrix12(size)
      real*8 matrix14(size)
      integer iseed
      integer readflag, readflag2
      real*8 sigpA(3), sigpB(3), sigth, sigph, XX0
      integer xdim2, ydim2, zdim2
      real*8 dm, dpt, drap
      real*8 mmin, mmax, ptmin, ptmax, rapmin, rapmax
      real*8 matrix51(size)
      common /haft/ fname, comment, xdim, ydim, zdim, dp, dth, dph,
     +              pmin, pmax, thmin, thmax, phmin, phmax,
     +              matrix2, matrix3, matrix8, matrix9,
     +              matrix10, matrix12, matrix14,
     +              iseed, sigpA, sigpB, sigth, sigph, XX0,
     +              fname2, comment2, xdim2, ydim2, zdim2,
     +              dm, dpt, drap,
     +              mmin, mmax, ptmin, ptmax, rapmin, rapmax,
     +              matrix51, readflag, readflag2
