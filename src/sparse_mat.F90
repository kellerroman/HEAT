module sparse_mat

!The Y12M documentation:
!
!                =================================
!                Documentation of subroutine Y12MA
!                =================================
!
!     1.  Purpose
!
!
!     Y12MA   solves   sparse  systems  of  linear  algebraic
!     equations by Gaussian elimination.  The subroutine is a
!     "black  box  subroutine"  designed to solve efficiently
!     problems which contain only one system  with  a  single
!     right  hand side. The number of the input parameters is
!     minimized. The user must assign values only to NN, NN1,
!     N,  Z,  A,  SNR,  RNR, IHA and B according to the rules
!     described in Section 2.4 (see below).  It is  extremely
!     easy  to  modify  the  subroutine  to  the cases: (a) a
!     sequence of systems with  the  same  matrix  is  to  be
!     solved (note that one system with many right hand sides
!     can be rewritten as a sequence of systems with the same
!     matrix),  (b)  a sequence of systems whose matrices are
!     different but of the same structure is to be solved and
!     (c)  a  sequence  of  systems whose matrices are of the
!     same structure and some of them are the same is  to  be
!     solved.   These  cases  are  defined as case (ii), case
!     (iii) and case (iv)  in  Section  1.1.  "Scope  of  the
!     Y12M".   The recommendations in Section 1.3.6 should be
!     followed in order to modify the  subroutine  Y12MA  for
!     the  above  cases.  If  a  sequence  of  systems  whose
!     matrices are of different structure (this case  appears
!     as  case  (v)  in  Section  1.1)  is  to be solved then
!     subroutine Y12MA should be called in  the  solution  of
!     each system in the sequence.
!
!     2.  Calling sequence and declaration of the parameters
!
!     The  subroutine  is  written in FORTRAN and it has been
!     extensively tested with the FOR and FTN compilers on  a
!     UNIVAC  1100/82  computer,  at  the  Regional Computing
!     Centre at the University of  Copenhagen  (RECKU).  Many
!     examples  have  been run on an IBM 3033 computer at the
!     Northern Europe University Computing Centre (NEUCC) and
!     on  a  CDC Cyber 173 computer at the Regional Computing
!     Centre  at  the  University  of  Aarhus  (RECAU).   Two
!     different  versions  are  available: a single precision
!     version named Y12MAE and  a  double  precision  version
!     named  Y12MAF.  The calls of these two versions and the
!     declarations of the parameters are as follows.
!
!     A). Single precision version: Y12MAE
!
!           SUBROUTINE Y12MAE(N, Z, A, SNR, NN, RNR, NN1,
!          1  PIVOT, HA, IHA, AFLAG, IFLAG, B, IFAIL)
!
!           REAL A(NN), PIVOT(N), B(N), AFLAG(8)
!           INTEGER SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
!           INTEGER N, Z, NN, NN1, IHA, IFAIL
!
!     B). Double precision version: Y12MAF
!
!           SUBROUTINE Y12MAF(N, Z, A, SNR, NN, RNR, NN1,
!          1  PIVOT, HA, IHA, AFLAG, IFLAG, B, IFAIL)
!
!           DOUBLE PRECISION A(NN), PIVOT(N), B(N), AFLAG(8)
!           INTEGER SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
!           INTEGER N, Z, NN, NN1, IHA, IFAIL
!
!     These two versions can be used on many other  computers
!     also. However some alterations may be needed and/or may
!     ensure greater efficiency of  the  performance  of  the
!     subroutine. For example, it will be much more efficient
!     to declare arrays SNR, RNR  and  (if  possible)  HA  as
!     INTEGER*2 arrays on some IBM installations.
!
!
!     3.  Method
!
!     The system Ax = b is solved  by  Gaussian  elimination.
!     Pivotal interchanges are used in an attempt to preserve
!     both the stability of the computations and the sparsity
!     of  the original matrix. In this way a decomposition LU
!     = PAQ is normally calculated. P and Q  are  permutation
!     matrices, L is a lower triangular matrix, U is an upper
!     triangular matrix. The right hand side vector b is also
!     modified  during  the  decomposition so that the vector
!     c=L**(-1)*P*b  is  available  after  the  decomposition
!     stage.  In  this  way  there  is  no  need to store the
!     non-zero elements of matrix L. Therefore these elements
!     are  not  stored.  This  leads  to  a  reduction of the
!     storage requirements (the length of arrays  A  and  SNR
!     can  be  decreased).  An  approximation to the solution
!     vector is found by solving U*tr(T)*x=c, where tr(T)  is
!     the  transpose.  The  subroutine calculates the rate at
!     which the  elements  of  the  matrix  grow  during  the
!     decomposition  (see  parameter  AFLAG(5) below) and the
!     minimal  in  absolute  value   pivotal   element   (see
!     parameter  AFLAG(8)  below).  These  two numbers can be
!     used to decide whether the  computed  approximation  is
!     acceptable   or  not.  Positive  values  of  a  special
!     parameter IFAIL indicate that the subroutine is  unable
!     to  solve  the problem. The error diagnostics, given in
!     Section 7, describe the probable cause.
!
!     4.  Parameters of the subroutine
!
!
!     N     - INTEGER On entry N must contain the  number  of
!             equations  in  the  system   Ax=b. Unchanged on
!             exit.
!
!     Z     - INTEGER. On entry Z must contain the number  of
!             non-zero  elements  in the coefficient matrix A
!             of the system Ax = b. Unchanged on exit.
!
!     A     - REAL (in the single precision  version  Y12MAE)
!             or  DOUBLE  PRECISION  (in the double precision
!             version Y12MAF) array of length NN (see below).
!             On  entry the first Z locations of array A must
!             contain   the   non-zero   elements   of    the
!             coefficient  matrix A of the system Ax = b. The
!             order  of  the   non-zero   elements   may   be
!             completely arbitrary. The content of array A is
!             modified by  subroutine  Y12MA.  On  successful
!             exit array A will contain the non-zero elements
!             of the upper triangular matrix U  (without  the
!             diagonal  elements  of  matrix  U  which can be
!             found in array PIVOT, see below).
!
!     SNR   - INTEGER array of  length  NN  (see  below).  On
!             entry  SNR(j),  j  =  1(1)Z,  must  contain the
!             column number of the non-zero element stored in
!             A(j).   The content of array SNR is modified by
!             subroutine Y12MA.  On successful exit array SNR
!             will contain the column numbers of the non-zero
!             elements  of  the  upper  triangular  matrix  U
!             (without  the  column  numbers  of the diagonal
!             elements of matrix U).
!
!     NN    - INTEGER. On entry NN must contain the length of
!             arrays   A  and  SNR.  Restriction:  NN.GE.2*Z.
!             Recommended value: 2*Z.LE.NN.LE.3*Z.  Unchanged
!             on exit.
!
!     RNR   - INTEGER  array  of  length  NN1 (see below). On
!             entry RNR(i), i = 1(1)Z, must contain  the  row
!             number  of the non-zero element stored in A(i).
!             The  content  of  array  RNR  is  modified   by
!             subroutine   Y12MA.   On  successful  exit  all
!             components of array RNR will normally be  zero.
!
!     NN1   - INTEGER. On entry NN1 must contain  the  length
!             of   the  array  RNR.   Restriction:  NN1.ge.Z.
!             Recommended      value:      2*Z.LE.NN1.LE.3*Z.
!             Unchanged on exit.
!
!     PIVOT - REAL  (in  the single precision version Y12MAE)
!             or DOUBLE PRECISION (in  the  double  precision
!             version  Y12MAF) array of length N. The content
!             of array PIVOT is modified by subroutine Y12MA.
!             On successful exit array PIVOT will contain the
!             pivotal  elements  (the  diagonal  elements  of
!             matrix  U). This means that a small element (or
!             small elements) in  array  PIVOT  on  exit  may
!             indicate    numerical    singularity   of   the
!             coefficient matrix A. Note that the smallest in
!             absolute  value  element in array PIVOT is also
!             stored in AFLAG(8), see below.
!
!     HA    - INTEGER two-dimensional array.  The  length  of
!             the  first  dimension  is  IHA (see below). The
!             length of the  second  dimension  is  11.  This
!             array  is  used  as  a work space by subroutine
!             Y12MA.
!
!     IHA   - INTEGER. On entry IHA must contain  the  length
!             of   the   first   dimension   of   array   HA.
!             Restriction:  IHA.GE.N. Unchanged on exit.
!
!     AFLAG - REAL (in the single precision  version  Y12MAE)
!             or  DOUBLE  PRECISION  (in the double precision
!             version Y12MAF) array of length 8. The  content
!             of array AFLAG is modified by subroutine Y12MA.
!             The content of the components of this array can
!             be described as follows.
!
!             AFLAG(1) - Stability factor. An element can  be
!                        chosen  as  pivotal  element only if
!                        this element is larger (in  absolute
!                        value)  than  the  absolute value of
!                        the  largest  element  in  its   row
!                        divided   by  AFLAG(1).   Subroutine
!                        Y12MA  sets  AFLAG(1)  =  16.   This
!                        value  has  been  found satisfactory
!                        for a wide range of test-examples.
!
!             AFLAG(2) - Drop-tolerance. An element, which in
!                        the   process  of  the  computations
!                        becomes smaller (in absolute  value)
!                        than  the drop-tolerance, is removed
!                        from array A (and its column and row
!                        numbers  are removed from arrays SNR
!                        and  RNR).   Subroutine  Y12MA  sets
!                        AFLAG(2)  =  1.0E-12. This value has
!                        been found satisfactory for  a  wide
!                        range   of  test-matrices.  By  this
!                        choice it is assumed that the matrix
!                        is not too badly scaled and that the
!                        magnitude of the elements is 1.
!
!             AFLAG(3) - The   subroutine   will   stop   the
!                        computations  when the growth factor
!                        (parameter  AFLAG(5),   see   below)
!                        becomes  larger  than  AFLAG(3). Our
!                        experiments show that if AFLAG(3)  >
!                        1.0E6  then the solution is normally
!                        quite wrong.   Therefore  subroutine
!                        Y12MA sets AFLAG(3) = 1.0E6.
!
!             AFLAG(4) - The   subroutine   will   stop   the
!                        computations when the absolute value
!                        of  a  current  pivotal  element  is
!                        smaller    than    AFLAG(4)*AFLAG(6)
!                        (where  the  absolute  value  of the
!                        largest (in absolute value)  element
!                        of  the original matrix is stored in
!                        AFLAG(6),    see     below).     Our
!                        experiments show that
!                        AFLAG(4)=1.0E-12  will normally be a
!                        good  choice.  Therefore  subroutine
!                        Y12MA sets  AFLAG(4)=1.0E-12.
!
!             AFLAG(5) - Growth factor. After each  stage  of
!                        the   elimination  subroutine  Y12MA
!                        sets  AFLAG(5)  =  AFLAG(7)/AFLAG(6)
!                        (the   description   of   parameters
!                        AFLAG(6)  and  AFLAG(7)   is   given
!                        below).  Large  values  of  AFLAG(5)
!                        indicate that an  appreciable  error
!                        in    the   computed   solution   is
!                        possible. In an  extreme  case  when
!                        AFLAG(5)  becomes larger than 1.0E6,
!                        see parameter  AFLAG(3),  subroutine
!                        Y12MA  will stop the computations in
!                        an attempt to prevent overflows.
!
!             AFLAG(6) - Subroutine Y12MA sets AFLAG(6) equal
!                        to  max(abs(a(i,j))), i = 1(1)N, j =
!                        1(1)N.
!
!             AFLAG(7) - Subroutine Y12MA sets AFLAG(7) equal
!                        to max(abs(a(i,j;s))),  1.lt.s.lt.k,
!                        after each step  k,   k=1(1)N-1,  of
!                        the elimination.
!
!             AFLAG(8) - Subroutine Y12MA sets AFLAG(8) equal
!                        to min(abs(a(i,i))), i= 1(1)N.  This
!                        means   that   the  minimal  pivotal
!                        element (in absolute value) will  be
!                        stored  in  AFLAG(8)  on  successful
!                        exit.   Small  values  of   AFLAG(8)
!                        indicate  numerical  singularity  of
!                        the original matrix.  We advise  the
!                        user to check this parameter on exit
!                        very carefully.
!
!     IFLAG - INTEGER array of length 10. The content of this
!             array is  modified  by  subroutine  Y12MA.  The
!             content  of the components of this array can be
!             described as follows.
!
!             IFLAG(1) - Subroutine  Y12MA uses IFLAG(1) as a
!                        work space.
!
!             IFLAG(2) - Subroutine Y12MA sets IFLAG(2) =  3.
!                        This means that at each stage of the
!                        Gaussian  elimination  (except   the
!                        last  one)  the  pivotal  search  is
!                        carried out in the 3 rows which have
!                        least  numbers of non-zero elements.
!
!             IFLAG(3) - Subroutine  Y12MA sets IFLAG(3) = 1.
!                        This means that the pivotal strategy
!                        for  general  matrices will be used.
!                        For  some  special   matrices   (for
!                        example positive definite) this will
!                        be  inefficient.   Subroutine  Y12MA
!                        can  easily  be  modified  for  such
!                        matrices;   only    the    statement
!                        IFLAG(3)  = 1 should be changed (for
!                        positive  definite   matrices   e.g.
!                        IFLAG(3)  =  2  can  be  used). More
!                        details  about  the  use   of   this
!                        parameter  with special matrices are
!                        given in Section 1.3.1
!
!             IFLAG(4) - Subroutine Y12MA sets IFLAG(4) =  0.
!                        This  is the best choice in the case
!                        where only one system Ax = b  is  to
!                        be   solved,  see  more  details  in
!                        Section 1.3.
!
!             IFLAG(5) - Subroutine Y12MA sets IFLAG(5) =  1.
!                        This   means   that   the   non-zero
!                        elements  of  the  lower  triangular
!                        matrix  L will be removed during the
!                        decomposition   stage,   see    more
!                        details in Section 1.3.3.
!
!             IFLAG(6) - On  successful exit IFLAG(6) will be
!                        equal to  the  number  of  "garbage"
!                        collections in the row ordered list.
!                        If IFLAG(6)  is  large  then  it  is
!                        better  to  choose a larger value of
!                        NN in the next calls  of  subroutine
!                        Y12MA  with  the  same  or a similar
!                        matrix  A.  This  will  lead  to   a
!                        reduction in the computing time.
!
!             IFLAG(7) - On  successful exit IFLAG(7) will be
!                        equal to  the  number  of  "garbage"
!                        collections  in  the  column ordered
!                        list.  If IFLAG(7) is large then  it
!                        is  better  to choose a larger value
!                        of  NN1  in  the   next   calls   of
!                        subroutine  Y12MA with the same or a
!                        similar matrix A. This will lead  to
!                        a reduction in the computing time.
!
!             IFLAG(8) - On  successful exit IFLAG(8) will be
!                        equal  to  the  maximal  number   of
!                        non-zero elements kept in array A at
!                        any   stage    of    the    Gaussian
!                        elimination.  If  IFLAG(8)  is  much
!                        smaller then NN (or  NN1)  then  the
!                        length  NN  (or  NN1)  can be chosen
!                        smaller  in  the   next   calls   of
!                        subroutine  Y12MA with the same or a
!                        similar matrix A. This will lead  to
!                        a reduction of the storage needed.
!
!             IFLAG(9) - This   parameter   is   ignored   by
!                        subroutine Y12MA.  It will  be  used
!                        in  some  of  the  other subroutines
!                        when IFLAG(4) = 1  is  specified  on
!                        entry.
!
!             IFLAG(10)- This   parameter   is   ignored   by
!                        subroutine Y12MA.  It will  be  used
!                        in  some  of  the  other subroutines
!                        when IFLAG(4) = 1  is  specified  on
!                        entry.
!
!     B     - REAL (in the single precision  version  Y12MAE)
!             or  DOUBLE  PRECISION  (in the double precision
!             version Y12MAF) array of length N. On entry the
!             right  hand  side vector b of the system Ax = b
!             must be stored in array B. The content of array
!             B   is   modified   by  subroutine  Y12MA.   On
!             successful exit the  computed  solution  vector
!             will be stored in array B.
!
!     IFAIL - Error  diagnostic  parameter.  The  content  of
!             parameter  IFAIL  is  modified  by   subroutine
!             Y12MA.  On exit IFAIL = 0 if the subroutine has
!             not detected  any  error.  Positive  values  of
!             IFAIL  on  exit  show  that some error has been
!             detected by the subroutine. Many of  the  error
!             diagnostics  are  common for all subroutines in
!             the package.  Therefore the  error  diagnostics
!             are listed in a separate section, Section 7, of
!             this book. We advise  the  user  to  check  the
!             value of this parameter on exit.
!
!     5.  Error diagnostics
!
!     Error diagnostics are given by positive values  of  the
!     parameter  IFAIL  (see  above).   We advise the user to
!     check carefully the value of this  parameter  on  exit.
!     The error messages are listed in Section 7.
!
!     6.  Auxiliary subroutines
!
!     Y12MA  calls  three other subroutines: Y12MB, Y12MC and
!     Y12MD.
!
!     7.  Timing
!
!     The time taken depends  on  the  order  of  the  matrix
!     (parameter  N),  the number of the non-zero elements in
!     the matrix (parameter Z), the magnitude of the non-zero
!     elements and their distribution in the matrix.
!
!     8.  Storage
!
!     There are no internally declared arrays.
!
!     9.  Accuracy
!
!     It  is  difficult  to  evaluate  the  accuracy  of  the
!     computed solution. Large values of parameter  AFLAG(5),
!     the   growth  factor,  indicate  unstable  computations
!     during  the  decomposition  stage.  Small   values   of
!     parameter  AFLAG(8)  can  be considered as a signal for
!     numerical singularity.  We  must  emphasize  here  that
!     normally much more reliable evaluations of the accuracy
!     achieved  can  be  found  by  the  use   of   iterative
!     refinement, i.e. by the use of subroutine Y12MF. By the
!     use of the latter subroutine the  computing  time  will
!     often be reduced too. However, the storage requirements
!     may be increased sometimes.
!
!     10.  Some remarks
!
!     Remark 1   Subroutine  Y12MA  may also be considered as
!                an example of how to initialize the first  4
!                components  of  array  AFLAG and the first 5
!                components of  array  IFLAG  when  only  one
!                system  with  a single right hand side is to
!                be  solved.  An   efficient   use   of   the
!                subroutines  Y12MB,  Y12MC and Y12MD for the
!                other cases discussed in Section 1.1 may  be
!                achieved  by  a single modification of Y12MA
!                according to  the  rules  given  in  Section
!                1.3.6   (a   change   of  only  one  or  two
!                statements in Y12MA is needed to obtain such
!                a modification).
!
!     Remark 2   The last 4 components of array AFLAG and the
!                last 5 components of array IFLAG can be used
!                as  output  parameters.   The  user  is  not
!                obliged to do this.  However,  we  recommend
!                the check of these parameters on exit.
!                =================================
!                Documentation of subroutine Y12MB
!                =================================
!
!     1.  Purpose
!
!     Y12MB prepares a system of linear  algebraic  equations
!     to  be  factorized (by subroutine Y12MC) and solved (by
!     subroutine Y12MD).  Each call of subroutine Y12MC  must
!     be  preceded  by a call of Y12MB. It is very convenient
!     to perform some operations on the non-zero elements  of
!     the  coefficient  matrix  between  calls  of  Y12MB and
!     Y12MC.  For  example,  when  iterative  refinement   is
!     carried  out  the information needed in the calculation
!     of the residual vector is copied in some  extra  arrays
!     between  the  calls of Y12MB and Y12MC. It is also very
!     easy to perform any kind of scaling after the  call  of
!     Y12MB.
!
!     2.  Calling sequence and declaration of the parameters
!
!     The  subroutine  is  written  in  FORTRAN  and has been
!     extensively tested with the FOR and FTN compilers on  a
!     UNIVAC  1100/82  computer  at  the  Regional  Computing
!     Centre at the University of  Copenhagen  (RECKU).  Many
!     examples  have  been run on an IBM 3033 computer at the
!     Northern Europe University Computing Centre (NEUCC) and
!     on  a  CDC Cyber 173 computer at the Regional Computing
!     Centre  at  the  University  of  Aarhus  (RECAU).   Two
!     different  versions  are available:  a single precision
!     version Y12MBE and a double precision  version  Y12MBF.
!     The  calls of these two versions and the declaration of
!     the parameters are as follows:
!
!     A) Single precision version Y12MBE
!
!           SUBROUTINE Y12MBE(N, Z, A, SNR, NN, RNR, NN1,
!          1  HA, IHA, AFLAG, IFLAG, IFAIL)
!
!           REAL A(NN), AFLAG(8)
!           INTEGER SNR(NN), RNR(NN1), HA(N,11), IFLAG(10)
!           INTEGER N, Z, IHA, IFAIL
!
!     B) Double precision version Y12MBF
!
!           SUBROUTINE Y12MBF(N, Z, A, SNR, NN, RNR, NN1,
!          1  HA, IHA, AFLAG, IFLAG, IFAIL)
!
!           DOUBLE PRECISION A(NN), AFLAG(8)
!           INTEGER SNR(NN), RNR(NN1), HA(N,11), IFLAG(10)
!           INTEGER N, Z, IHA, IFAIL
!
!     These  two versions can be used on many other computers
!     also. However some alterations may be needed and/or may
!     ensure  greater  efficiency  of  the performance of the
!     subroutine. For example, it will be much more efficient
!     to  declare  arrays  SNR,  RNR  and (if possible) HA as
!     INTEGER*2 arrays on some IBM installations.
!
!
!     3.  Method
!
!     The  non-zero  elements of matrix A are ordered by rows
!     and stored in the first Z positions  of  array  A.  The
!     order   of  the  non-zero  elements  within  a  row  is
!     arbitrary. The column numbers of the non-zero  elements
!     are  stored  in  the  first Z positions of array SNR so
!     that if  A(J)=a(i,j),  J=1(1)Z, then  SNR(J)=j. The row
!     numbers  of  the  non-zero  elements  are stored in the
!     first Z positions of array RNR so that the row  numbers
!     of  the non-zero elements in the first column of matrix
!     A are located before the row numbers  of  the  non-zero
!     elements  in  the  second  column  of matrix A, the row
!     numbers of the non-zero elements in the  second  column
!     of  matrix  A are located before the row numbers of the
!     non-zero elements in the third column of matrix  A  and
!     so  on. Some additional information, e.g. about the row
!     starts, row ends, column starts  and  column  ends,  is
!     stored  in  the  work array HA. This storage scheme has
!     been proposed by Gustavson.
!
!     4.  Parameters of the subroutine
!
!
!     N     - INTEGER. On entry N must contain the number  of
!             equations  in  the  system   Ax=b. Unchanged on
!             exit.
!
!     Z     - INTEGER. On entry Z must contain the number  of
!             non-zero  elements  in the coefficient matrix A
!             of the system Ax = b. Unchanged on exit.
!
!     A     - REAL (in the single precision  version  Y12MBE)
!             or  DOUBLE  PRECISION  (in the double precision
!             version Y12MBF) array of length NN (see below).
!             On  entry the first Z locations of array A must
!             contain   the   non-zero   elements   of    the
!             coefficient  matrix A of the system Ax = b. The
!             order  of  the   non-zero   elements   may   be
!             completely arbitrary. The content of array A is
!             modified by  subroutine  Y12MB.  On  successful
!             exit  the  first  Z  positions  of array A will
!             contain the non-zero elements of the  matrix  A
!             ordered by rows.
!
!     SNR   - INTEGER  array  of  length  NN  (see below). On
!             entry SNR(j),  j  =  1(1)Z,  must  contain  the
!             column number of the non-zero element stored in
!             A(j).  The content of array SNR is modified  by
!             subroutine  Y12MB.   On successful exit SNR(j),
!             j=1(1)Z, will contain the column number of  the
!             non-zero element stored in A(j).
!
!     NN    - INTEGER. On entry NN must contain the length of
!             arrays  A  and  SNR.  Restriction:   NN.ge.2*Z.
!             Recommended  value: 2*Z.le.NN.le.3*Z.  See also
!             the description of NN in Y12MC.   Unchanged  on
!             exit.
!
!     RNR   - INTEGER  array  of  length  NN1 (see below). On
!             entry RNR(i), i = 1(1)Z, must contain  the  row
!             number  of the non-zero element stored in A(i).
!             The  content  of  array  RNR  is  modified   by
!             subroutine  Y12MB.  On  successful exit the row
!             numbers of the non-zero elements  of  matrix  A
!             are  stored  in  the first Z positions of array
!             RNR, so that the row numbers  of  the  non-zero
!             elements  in  the  first column of matrix A are
!             before the row numbers of the non-zero elements
!             in  the  second  column  of  matrix  A, the row
!             numbers of the non-zero elements in the  second
!             column  of  matrix A are before the row numbers
!             of the non-zero elements in the third column of
!             matrix A and so on. This means that on exit the
!             row number of the non-zero  element  stored  in
!             A(i),  i  =  1(1)Z, is not stored in RNR(i), in
!             general.
!
!     NN1   - INTEGER. On entry NN1 must contain  the  length
!             of   the  array  RNR.   Restriction:  NN1.ge.Z.
!             Recommended value: 2*Z.le.NN1.le.3*Z.  See also
!             the  desription  of NN1 in Y12MC.  Unchanged on
!             exit.
!
!     HA    - INTEGER two-dimensional array.  The  length  of
!             the  first  dimension of HA is IHA (see below).
!             The length of the second dimension of HA is 11.
!             The contents of some elements of this array are
!             modified by subroutine Y12MB, the  contents  of
!             the others are ignored. On successful exit:
!
!             (1) HA(i,1) contains the position  in  array  A
!                         where the first element of row i, i
!                         = 1(1)N, is stored.
!             (2) HA(i,3) contains  the  position  in array A
!                         where the last element of row i,  i
!                         =  1(1)N,  is  stored (all non-zero
!                         elements of  row   i   are  located
!                         between    HA(i,1)    and   HA(i,3)
!                         compactly; the number  of  non-zero
!                         elements  in  row  i  is  HA(i,3) -
!                         HA(i,1) + 1).
!             (3) HA(i,4) contains  the position in array RNR
!                         where the row number of  the  first
!                         non-zero  element  of  column  j is
!                         stored (j = 1(1)N).
!             (4) HA(j,6) contains  the position in array RNR
!                         where the row number  of  the  last
!                         element  in column j, j = 1(1)N, is
!                         stored  (all  row  numbers  of  the
!                         non-zero  elements  of column j are
!                         located between HA(j,4) and HA(j,6)
!                         compactly,  the  number of non-zero
!                         elements    in    column    j    is
!                         HA(j,6)-HA(j,4)+1).
!
!             (5)         Some   information  needed  in  the
!                         pivotal    search    during     the
!                         decomposition  is stored in columns
!                         7,8 and 11 of array HA.
!
!             (6)         The other columns of HA are  either
!                         ignored  by Y12MB or used as a work
!                         space.
!
!     IHA   - INTEGER.  On  entry IHA must contain the length
!             of   the   first   dimension   of   array   HA.
!             Restriction:  IHA.ge.N. Unchanged on exit.
!
!     AFLAG - REAL  (in  the single precision version Y12MBE)
!             or DOUBLE PRECISION (in  the  double  precision
!             version Y12MBF) array of length 8. The contents
!             of  all  locations  of  array   AFLAG,   except
!             AFLAG(6), are ignored by subroutine Y12MB.  The
!             content of AFLAG(6) is modified  by  Y12MB.  On
!             successful exit AFLAG(6) contains
!             max(abs(i,j)).
!
!     IFLAG - INTEGER array of length 10. The user should set
!             IFLAG(1).ge.0  before  the call of package Y12M
!             (i.e. before the first call of a subroutine  of
!             this  package).  IFLAG(1)  is used in the error
!             checks by the subroutines  and  should  not  be
!             altered  by the user between any two successive
!             calls of subroutines of the  package.  IFLAG(1)
!             will  be  equal to  -1  on successful exit from
!             subroutine Y12MB. Thus, the content of IFLAG(1)
!             is  modified  by  subroutine  Y12MB.  On  entry
!             IFLAG(4) must contain 0, 1 or 2. IFLAG(4)  =  0
!             is  the  best choice when only one system is to
!             be solved, when the first system of a  sequence
!             of systems with the same matrix is to be solved
!             and when any system of a  sequence  of  systems
!             whose matrices are of different structure is to
!             be solved. IFLAG(4) = 1 is the best choice when
!             the  first system in a sequence of systems with
!             the same structure is to be solved. IFLAG(4)  =
!             2  is the best choice when any system after the
!             first  one  in  a  sequence  of  systems  whose
!             matrices  are  of  the  same structure is to be
!             solved. The content of IFLAG(4) is unchanged by
!             subroutine  Y12MB. The other locations of IFLAG
!             are either ignored by Y12MB or used as  a  work
!             space.
!
!     IFAIL - Error  diagnostic  parameter.  The  content  of
!             parameter  IFAIL  is  modified  by   subroutine
!             Y12MB.  On exit IFAIL = 0 if the subroutine has
!             not detected  any  error.  Positive  values  of
!             IFAIL  on  exit  show  that some error has been
!             detected by the subroutine. Many of  the  error
!             diagnostics  are  common for all subroutines in
!             the package.  Therefore the  error  diagnostics
!             are listed in a separate section, Section 7, of
!             this book. We advise  the  user  to  check  the
!             value of this parameter on exit.
!
!     5.  Error diagnostics
!
!     Error diagnostics are given by positive values  of  the
!     parameter  IFAIL  (see  above).  We  advise the user to
!     check carefully the value of this  parameter  on  exit.
!     The error messages are listed in Section 7.
!
!     6.  Auxiliary subroutines
!
!     None.
!
!     7.  Timing
!
!     The  time  taken  depends  on  the  number  of non-zero
!     elements (parameter Z) in matrix A.
!
!     8.  Storage
!
!     There are no internally declared arrays.
!
!     9.  Some remarks
!
!     Remark 1 The use of subroutine Y12MB is followed by the
!              use of Y12MC and Y12MD. Subroutines Y12MA  and
!              Y12MF  can be considered as examples of how to
!              use Y12MB, Y12MC and Y12MD.
!
!     Remark 2 The contents of N, Z, A, SNR,  NN,  RNR,  NN1,
!              columns  1,  3,  4, 6, 7, 8 and 11 of HA, IHA,
!              AFLAG(6), IFLAG(1), IFLAG(4) and IFAIL  should
!              not  be  altered  between  calls  of Y12MB and
!              Y12MC.
!
!     Remark 3 If IFAIL > 0 on exit  of  Y12MB  there  is  no
!              sense  in  calling  Y12MC and Y12MD. Therefore
!              the computations should  be  stopped  in  this
!              case  and  one  should  investigate (using the
!              error diagnostics  given  in  Section  7)  why
!              Y12MB  has assigned a positive value to IFAIL.
!                =================================
!                Documentation of subroutine Y12MC
!                =================================
!
!     1.  Purpose
!
!     Y12MC   decomposes  a  matrix  A  into  two  triangular
!     matrices L and U. Each call of subroutine Y12MC  should
!     be  preceded  by  a  call of Y12MB. Subroutine Y12MD is
!     normally called one or several times after the call  of
!     Y12MC.
!
!     2.  Calling sequence and declaration of the parameters
!
!     The  subroutine  is  written  in  FORTRAN  and has been
!     extensively tested with the FOR and  FTN  compilers  on
!     the  UNIVAC  1100/82 computer at the Regional Computing
!     Centre at the University of  Copenhagen  (RECKU).  Many
!     examples  have  been run on an IBM 3033 computer at the
!     Northern Europe University Computing Centre (NEUCC) and
!     on  a  CDC Cyber 173 computer at the Regional Computing
!     Centre  at  the  University  of  Aarhus  (RECAU).   Two
!     different  versions  are  available: a single precision
!     version named Y12MCE and  a  double  precision  version
!     named  Y12MCF.  The calls of these two versions and the
!     declaration of the parameters are as follows.
!
!     A). Single precision version: Y12MCE
!
!           SUBROUTINE Y12MCE(N, Z, A, SNR, NN, RNR, NN1,
!          1  PIVOT, B, HA, IHA, AFLAG, IFLAG, IFAIL)
!
!           REAL A(NN), PIVOT(N), B(N), AFLAG(8)
!           INTEGER SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
!           INTEGER N, Z, NN, NN1, IHA, IFAIL
!
!     B). Double precision version: Y12MCF
!
!           SUBROUTINE Y12MCF(N, Z, A, SNR, NN, RNR, NN1,
!          1  PIVOT, B, HA, IHA, AFLAG, IFLAG, IFAIL)
!
!           DOUBLE PRECISION A(NN), PIVOT(N), B(N), AFLAG(8)
!           INTEGER SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
!           INTEGER N, Z, NN, NN1, IHA, IFAIL
!
!     These  two versions can be used on many other computers
!     also. However some alterations may be needed and/or may
!     ensure  greater  efficiency  of  the performance of the
!     subroutine. For example, it will be much more efficient
!     to  declare  arrays  SNR,  RNR  and (if possible) HA as
!     INTEGER*2 arrays on some IBM installations.
!
!
!     3.  Method
!
!     Gaussian  elimination  is  used in the decomposition of
!     matrix A. Pivotal interchanges are  implemented  as  an
!     attempt   to   preserve   both  the  stability  of  the
!     computations  and  the  sparsity  of  the   coefficient
!     matrix. Two triangular matrices L and U are computed so
!     that  LU=PAQ  (where P and Q are permutation matrices).
!     The  right  hand side vector b is also modified so that
!     vector c = L**(-1)*P*b is computed on successful  exit.
!     In  this  way  matrix  L is sometimes not needed in the
!     further computations and may be removed if this is  so.
!
!     4.  Parameters of the subroutine
!
!
!     N     - INTEGER. On entry N must contain the number  of
!             equations  in  the  system   Ax=b. Unchanged on
!             exit.
!
!     Z     - INTEGER. On entry Z must contain the number  of
!             non-zero  elements  in the coefficient matrix A
!             of the system  Ax=b. Unchanged on exit.
!
!     A     - REAL (in the single precision  version  Y12MCE)
!             or  DOUBLE  PRECISION  (in the double precision
!             version Y12MCF) array of length NN (see below).
!             On  entry the first Z locations of array A must
!             contain   the   non-zero   elements   of    the
!             coefficient  matrix  A  of  the  system  Ax = b
!             ordered by  rows  (in  the  preceding  call  of
!             Y12MB).  The  content of array A is modified by
!             subroutine Y12MC.  On successful exit  array  A
!             will contain the non-zero elements of the upper
!             triangular  matrix  U  (without  the   diagonal
!             elements  of  matrix  U  which  can be found in
!             array PIVOT, see below) and sometimes also  the
!             non-zero   elements  of  the  lower  triangular
!             matrix  L (without the diagonal elements  which
!             are not stored).
!
!     SNR   - INTEGER  array  of  length  NN  (see below). On
!             entry SNR(j),  j  =  1(1)Z,  must  contain  the
!             column number of the non-zero element stored in
!             A(j).   This  is  accomplished  by   subroutine
!             Y12MB.  The content of array SNR is modified by
!             subroutine Y12MC.  On successful exit array SNR
!             will contain the column numbers of the non-zero
!             elements  of  the  upper  triangular  matrix  U
!             (without  the  column  numbers  of the diagonal
!             elements of matrix U) and  sometimes  also  the
!             column  numbers of the non-zero elements of the
!             lower triangular matrix L (without  the  column
!             numbers of the diagonal elements of L).
!
!     NN    - INTEGER. On entry NN must contain the length of
!             array  A  and  SNR.   Restriction:   NN.ge.2*Z.
!             Recommended  value:  2*Z .le.NN.le.  3*Z if the
!             non-zero elements of matrix L will  be  removed
!             by  subroutine  Y12MC  (i.e. if Y12MC is called
!             with IFLAG(5) = 1) or 3*Z.le.  NN.le.5*Z if the
!             non-zero  elements  of matrix L will be kept by
!             subroutine Y12MC (i.e. if Y12MC is called  with
!             IFLAG(5) = 2).  Unchanged on exit.
!
!     RNR   - INTEGER  array  of  length  NN1 (see below). On
!             entry the first Z locations of array  RNR  must
!             contain   the   row  numbers  of  the  non-zero
!             elements of the coefficient  matrix  A  of  the
!             system  Ax  =  b so that the row numbers of the
!             non-zero elements in the first column of matrix
!             A  are  located  before  the row numbers of the
!             non-zero  elements  in  the  second  column  of
!             matrix  A,  the  row  numbers  of  the non-zero
!             elements in the second column of matrix  A  are
!             located  before the row numbers of the non-zero
!             elements in the third column of matrix A and so
!             on  (this  ordereing  is prepared by subroutine
!             Y12MB). The content of array RNR is modified by
!             subroutine   Y12MC.   On  successful  exit  all
!             components of array RNR will normally be  zero.
!
!     NN1   - INTEGER. On entry NN1 must contain  the  length
!             of   the  array  RNR.   Restriction:  NN1.ge.Z.
!             Recommended     value:       2*Z.le.NN1.le.3*Z.
!             Unchanged on exit.
!
!     PIVOT - REAL  (in  the single precision version Y12MCE)
!             or DOUBLE PRECISION (in  the  double  precision
!             version  Y12MCF) array of length N. The content
!             of array PIVOT is modified by subroutine Y12MC.
!             On successful exit array PIVOT will contain the
!             pivotal elements (  the  diagonal  elements  of
!             matrix  U). This means that a small element (or
!             small elements) in  array  PIVOT  on  exit  may
!             indicate    numerical    singularity   of   the
!             coefficient matrix A. Note that the smallest in
!             absolute value element in array PIVOT is stored
!             in AFLAG(8), see below.
!
!     B -     REAL (in the single precision  version  Y12MCE)
!             or  DOUBLE  PRECISION  (in the double precision
!             version Y12MCF) array of length  N.  The  right
!             hand side vector b of the system Ax = b must be
!             stored in B on entry. The content of array B is
!             modified  by  subroutine  Y12MC.  On successful
!             exit array B will  contain  the  components  of
!             vector c = L**(-1)*P*b.
!
!     HA    - INTEGER  two-dimensional  array.  The length of
!             the first dimension of HA is IHA  (see  below).
!             The length of the second dimension of HA is 11.
!             The content of columns 1, 3, 4, 6, 7, 8 and  11
!             is  prepared by Y12MB (and thus they should not
!             be altered between calls of Y12MB  and  Y12MC).
!             The   content   of  array  HA  is  modified  by
!             subroutine   Y12MC.   The   non-zero   elements
!             (without  the  diagonal  elements)  of row i in
!             matrix L  (when  this  matrix  is  stored)  are
!             between HA(i,1) and HA(i,2) -1. If the non-zero
!             elements  of  matrix  L  are  not  stored  then
!             HA(i,1)=HA(i,2), i=1(1)N. The non-zero elements
!             of row i in  matrix  U  (without  the  diagonal
!             elements)   are   stored  between  HA(i,2)  and
!             HA(i,3).  Information concerning  the  row  and
!             column  interchanges  is  stored in HA(i,7) and
!             HA(i,8),   i=1(1)N-1.   Information  about  the
!             largest  number  of  non-zero elements found in
!             row i /  column  j  during  any  stage  of  the
!             elimination  is  stored (when IFLAG(4)=1  only)
!             in  HA(i,9)/HA(j,10).  The  other   information
!             stored  in  array HA is not used in the further
!             computations.
!
!     IHA   - INTEGER. On entry IHA must contain  the  length
!             of   the   first   dimension   of   array   HA.
!             Restriction: IHA
!            .ge.N. Unchanged on exit.
!
!     AFLAG - REAL  (in  the single precision version Y12MCE)
!             or DOUBLE PRECISION (in  the  double  precision
!             version  Y12MCF) array of length 8. The content
!             of the array can be described as follows:
!
!             AFLAG(1) -  Stability factor. An element can be
!                         chosen as pivotal element  only  if
!                         this element is larger (in absolute
!                         value) than the absolute  value  of
!                         the  largest  element  in  its  row
!                         divided  by  AFLAG(1).   On   entry
!                         AFLAG(1)   should  contain  a  real
!                         number larger than 1.0. If this  is
!                         not  so  then  the  subroutine sets
!                         AFLAG(1)  =   1.0005.   Recommended
!                         values:   AFLAG(1) ranging from 4.0
!                         to 16.0.  Unchanged on  exit  (when
!                         correctly initialized).
!
!             AFLAG(2)  - Drop-tolerance. An element which in
!                         the  process  of  the  computations
!                         becomes smaller (in absolute value)
!                         than the drop-tolerance is  removed
!                         from  array  A  (and  its  row  and
!                         column  numbers  are  removed  from
!                         arrays   RNR  and  SNR).  On  entry
!                         AFLAG(2)   should   contain    some
!                         positive   small  number  or  zero.
!                         Recommended   value   AFLAG(2)    =
!                         1.0E-12  when  matrix  A is not too
!                         badly scaled and the  magnitude  of
!                         the   elements   is   1,  otherwise
!                         smaller values of  AFLAG(2)  should
!                         be used.  Unchanged on exit.
!
!              AFLAG(3) - The   subroutine   will   stop  the
!                         computation when the growth  factor
!                         (parameter   AFLAG(5),  see  below)
!                         becomes larger  than  AFLAG(3).  On
!                         entry  AFLAG(3)  should  contain  a
!                         large    positive    number.     If
!                         AFLAG(3)<1.0e5  then the subroutine
!                         sets  AFLAG(3)=1.0e5.   Recommended
!                         value   AFLAG(3)=1.0E6.   Unchanged
!                         on     exit     (when     correctly
!                         initialized).
!
!              AFLAG(4) - The   subroutine   will   stop  the
!                         computation when the absolute value
!                         of  a  current  pivotal  element is
!                         smaller   than    AFLAG(4)*AFLAG(6)
!                         (parameter  AFLAG(6)  is  described
!                         below).  On  entry  AFLAG(4)   must
!                         contain    a   small   non-negative
!                         number.  If  AFLAG(4)<0   on  entry
!                         then     the     subroutine    sets
!                         AFLAG(4)=-AFLAG(4).     Recommended
!                         value  AFLAG(4)=1.0E-12.  Unchanged
!                         on     exit     (when     correctly
!                         initialized).
!
!             AFLAG(5)  - Growth   factor.   The  content  of
!                         parameter AFLAG(5) is  modified  by
!                         subroutine  Y12MC. After each stage
!                         of   the    Gaussian    elimination
!                         subroutine  Y12MC  sets  AFLAG(5) =
!                         AFLAG(7)/AFLAG(6)       (parameters
!                         AFLAG(6) and AFLAG(7) are described
!                         below). On  exit  large  values  of
!                         parameters  AFLAG(5)  indicate that
!                         an   appreciable   error   in   the
!                         computed  solution  is possible. In
!                         an    extreme    case    ,    where
!                         AFLAG(5)>AFLAG(3)   the  subroutine
!                         will terminate the computations  in
!                         an attempt to prevent overflow (and
!                         IFAIL will be set equal to 4).
!
!             AFLAG(6)  - On entry AFLAG(6) is equal  to  the
!                         largest  element in the coefficient
!                         matrix A of the system Ax = b  (set
!                         by subroutine Y12MB).  Unchanged on
!                         exit.
!
!             AFLAG(7)  - On exit the  largest  (in  absolute
!                         value)  element  found  during  any
!                         stage of the  elimination  will  be
!                         stored  in AFLAG(7). The content of
!                         parameter AFLAG(7) is  modified  by
!                         subroutine Y12MC.
!
!             AFLAG(8)  - On  entry  the minimal (in absolute
!                         value)  pivotal  element  will   be
!                         stored in AFLAG(8). Small values of
!                         AFLAG(8)     indicate     numerical
!                         singularity   of   the  coefficient
!                         matrix A. We  advise  the  user  to
!                         check  this  parameter on exit from
!                         the calculation very carefully. The
!                         content  of  parameter  AFLAG(8) is
!                         modified by the subroutine Y12MC.
!
!     IFLAG - INTEGER array of length 10. The content of this
!             array is  modified  by  subroutine  Y12MC.  The
!             contents of the components of this array can be
!             described as follows.
!
!             IFLAG(1) - This parameter is used in connection
!                        with the error diagnostics. The user
!                        should  set IFLAG(1).ge.0 before the
!                        call of package  Y12M  (i.e.  before
!                        the  first  call  of a subroutine of
!                        this package). IFLAG(1) is  used  in
!                        the  error checks by the subroutines
!                        and should not  be  altered  by  the
!                        user   between  any  two  successive
!                        calls of subroutines of the package.
!                        IFLAG(1)  will  be  equal  to  -2 on
!                        successful  exit   from   subroutine
!                        Y12MC. Thus, the content of IFLAG(1)
!                        is modified by subroutine Y12MC.
!
!             IFLAG(2) - On entry IFLAG(2) must contain  some
!                        positive  integer smaller than N. We
!                        recommend IFLAG(2).le.3. If IFLAG(3)
!                        =  0  then this parameter is ignored
!                        by     subroutine     Y12MC.      If
!                        IFLAG((3).ge.0   then   the  pivotal
!                        search   at   any   stage   of   the
!                        elimination (except possibly some of
!                        the last stages) is carried  out  in
!                        the  IFLAG(2)  rows which have least
!                        numbers   of   non-zero    elements.
!                        Unchanged on exit.
!
!             IFLAG(3) - On  entry IFLAG(3) must contain 0, 1
!                        or 2.  For  general  pivotal  search
!                        IFLAG(3)  should  be set equal to 1.
!                        If IFLAG(3) = 2 then  only  diagonal
!                        elements of the coefficient matrix A
!                        will   be   selected   as    pivotal
!                        elements.  If  IFLAG(3) = 0 then the
!                        Gaussian elimination will be carried
!                        out     without     any    pivoting.
!                        IFLAG(3)=0 or IFLAG(3)=2  (i.e.  one
!                        of the special pivotal strategies is
!                        to be applied) should be  used  very
!                        carefully    because    the    error
!                        diagnostics algorithm may not  trace
!                        all   errors   in  this  case.   For
!                        example, if the user attempts to use
!                        IFLAG(3)=0   for   matrix   A  which
!                        contains zero elements on  the  main
!                        diagonal, then the run will often be
!                        stopped because a division  by  zero
!                        occurs.  Unchanged on exit.
!
!             IFLAG(4) - On  entry IFLAG(4) must contain 0, 1
!                        or 2.  IFLAG(4)  =  0  is  the  best
!                        choice  when  (i) only one system is
!                        to be solved, (ii) the first  system
!                        of  a  sequence  of systems with the
!                        same matrix (Ax = b1, Ax = b2,  ...,
!                        Ax = bp) is to be solved, (iii) when
!                        any system in a sequence of  systems
!                        whose   matrices  are  of  different
!                        structure is to be solved.  IFLAG(4)
!                        =  1  is  the  best  choice when the
!                        first  system  of  a   sequence   of
!                        systems  whose  matrices  are of the
!                        same structure is to be  solved.  In
!                        this case IFLAG(4) = 2 is to be used
!                        in the solution of any system  after
!                        the first one.  Unchanged on exit.
!
!             IFLAG(5) - On  entry IFLAG(5) must contain 1 or
!                        2.   If  IFLAG(5)  =  1   then   the
!                        non-zero  elements  of matrix L will
!                        be removed. If IFLAG(5) = 2 then the
!                        non-zero  elements  of matrix L will
!                        be stored.  Unchanged on exit.
!
!             IFLAG(6) - On successful exit IFLAG(6) will  be
!                        equal  to  the  number  of "garbage"
!                        collections in the row ordered list.
!                        If  IFLAG(6)  is  large  then  it is
!                        better to choose a larger  value  of
!                        NN  with  next  calls  of subroutine
!                        Y12MC with the  same  or  a  similar
!                        matrix   A.  This  will  lead  to  a
!                        reduction  in  the  computing  time.
!                        The  content of IFLAG(6) is modified
!                        by the subroutine Y12MC.
!
!             IFLAG(7) - On successful exit IFLAG(7) will  be
!                        equal  to  the  number  of "garbage"
!                        collections in  the  column  ordered
!                        list.   If IFLAG(7) is large then it
!                        is better to choose a  larger  value
!                        of   NN1   in   the  next  calls  of
!                        subroutine Y12MC with the same or  a
!                        similar  matrix A. This will lead to
!                        a reduction in the  computing  time.
!                        The  content of IFLAG(7) is modified
!                        by subroutine Y12MC.
!
!             IFLAG(8) - On successful exit IFLAG(8) will  be
!                        equal   to  the  maximal  number  of
!                        non-zero elements kept in array A at
!                        any    stage    of    the   Gaussian
!                        elimination.  If  IFLAG(8)  is  much
!                        smaller  then  NN  (or NN1) then the
!                        length NN (or  NN1)  can  be  chosen
!                        smaller  in next calls of subroutine
!                        Y12MC with the  same  or  a  similar
!                        matrix   A.  This  will  lead  to  a
!                        reduction of the storage needed. The
!                        content  of  IFLAG(8) is modified by
!                        subroutine Y12MC.
!
!             IFLAG(9) - The content of IFLAG(3) is  modified
!                        by  subroutine Y12MC when IFLAG(4) =
!                        1  and   ignored   otherwise.    The
!                        minimal  length  NN1 such that Y12MC
!                        can solve systems whose matrices are
!                        of   the   same   structure  without
!                        "garbage" collections in the  column
!                        ordered   list   and   "movings"  of
!                        columns at the  end  of  the  column
!                        ordered  list  is stored in IFLAG(9)
!                        after  the  solution  of  the  first
!                        system   in   the   sequence   (with
!                        IFLAG(4) = 1).
!
!             IFLAG(10)- The content of IFLAG(10) is modified
!                        by   the   subroutine   Y12MC   when
!                        IFLAG(4) = 1 and ignored  otherwise.
!                        The  minimal  length  NN  such  that
!                        subroutine Y12MC can  solve  systems
!                        whose   matrices  are  of  the  same
!                        structure     without      "garbage"
!                        collections  in the row ordered list
!                        and "movings" of rows to the end  of
!                        the  row  ordered  list is stored in
!                        IFLAG(10) after the solution of  the
!                        first  system  in the sequence (with
!                        IFLAG(4) = 1).
!
!     IFAIL - Error  diagnostics  parameter.  The  content of
!             parameter  IFAIL  is  modified  by   subroutine
!             Y12MC.  On exit IFAIL = 0 if the subroutine has
!             not detected  any  error.  Positive  values  of
!             IFAIL  on  exit  show  that some error has been
!             detected by the subroutine. Many of  the  error
!             diagnostics  are  common for all subroutines in
!             the package.  Therefore the  error  diagnostics
!             are listed in a separate section, Section 7, of
!             this book. We advise  the  user  to  check  the
!             value of this parameter on exit.
!
!     5.  Error diagnostics
!
!     Error diagnostics are given by positive values  of  the
!     parameter  IFAIL  (see  above).   We advise the user to
!     check carefully the value of this  parameter  on  exit.
!     The error messages are listed in Section 7.
!
!     6.  Auxiliary subroutines
!
!     None.
!
!     7.  Timing
!
!     The  time  taken  depends  on  the  order of the matrix
!     (parameter N), the number of the non-zero  elements  in
!     the matrix (parameter Z), the magnitude of the non-zero
!     elements and their distribution in the matrix.
!
!     8.  Storage
!
!     There are no internally declared arrays.
!
!     9.  Accuracy
!
!     It  is  difficult  to  evaluate  the  accuracy  of  the
!     computed  solution. Large values of parameter AFLAG(5),
!     the  growth  factor,  indicate  unstable   computations
!     during   the   decomposition  stage.  Small  values  of
!     parameter AFLAG(8) can be considered as  a  signal  for
!     numerical  singularity.  We  must  emphasize  here that
!     normally much more reliable evaluations of the accuracy
!     achieved   can   be  found  by  the  use  of  iterative
!     refinement, i.e. by the use of subroutine Y12MF. By the
!     use  of  the  latter subroutine the computing time will
!     often be reduced too.
!
!     10.  Some remarks
!
!     Remark 1   The values on entry of Y12MC are the same as
!                the  values  on  exit  of  Y12MB  for   many
!                parameters.  Therefore  the content of N, Z,
!                A, SNR, NN, RNR, NN1, columns 1, 3, 4, 6, 7,
!                8 and 11 of HA, AFLAG(6), IFLAG(1), IFLAG(4)
!                and IFAIL  should  not  be  changed  between
!                calls of Y12MB and Y12MC.
!
!     Remark 2   The  call  of  Y12MC is normally followed by
!                one or several calls of Y12MD.  The  content
!                of N, A, SNR, NN, B, PIVOT, columns 1, 2, 3,
!                7  and  8  of  HA,  IHA,  AFLAG,   IFLAG(1),
!                IFLAG(3),  IFLAG(4)  and IFAIL should not be
!                altered between calls of Y12MC and Y12MD.
!
!     Remark 3   If IFAIL > 0 on exit from Y12MC there is  no
!                justification  for  calling Y12MD. Therefore
!                the computations should  be  terminated  and
!                the user should investigate (using the error
!                diagnostics given in Section  7)  why  Y12MC
!                has assigned a positive value to IFAIL.
!                =================================
!                Documentation of subroutine Y12MD
!                =================================
!
!     1.  Purpose
!
!     Y12MD   solves   sparse  systems  of  linear  algebraic
!     equations.
!
!     2.  Calling sequence and declaration of the parameters
!
!     The subroutine is written in FORTRAN and  it  has  been
!     extensively  tested  with  the FOR and FTN compilers on
!     the UNIVAC 1100/82 computer at the  Regional  Computing
!     Centre  at  the  University of Copenhagen (RECKU). Many
!     examples have been run on an IBM 3033 computer  at  the
!     Northern Europe University Computing Centre (NEUCC) and
!     on a CDC Cyber 173 computer at the  Regional  Computing
!     Centre  at  the  University  of  Aarhus  (RECAU).   Two
!     different versions are available:  a  single  precision
!     version  named  Y12MDE  and  a double precision version
!     named Y12MDF. The calls of these two versions  and  the
!     declarations of the parameters are as follows.
!
!     A). Single precision version: Y12MDE
!
!           SUBROUTINE Y12MDE(N, A, NN, B, PIVOT, SNR,
!          1                  HA, IHA,, IFLAG, IFAIL)
!           REAL A(NN), PIVOT(N), B(N)
!           INTEGER SNR(NN), HA(IHA,11), IFLAG(10)
!           INTEGER N, NN, IHA, IFAIL
!
!     B). Double precision version: Y12MDF
!
!           SUBROUTINE Y12MDF(N, A, NN, B, PIVOT, SNR,
!          1                  HA, IHA,, IFLAG, IFAIL)
!           DOUBLE PRECISION A(NN), PIVOT(N), B(N)
!           INTEGER SNR(NN), HA(IHA,11), IFLAG(10)
!           INTEGER N, NN, IHA, IFAIL
!
!     These two versions can be used on many others computers
!     also. However some alterations may be needed and/or may
!     ensure  greater  efficiency  of  the performance of the
!     subroutine. For example, it will be much more efficient
!     to declare arrays SNR and (if possible) HA as INTEGER*2
!     arrays on some IBM installations.
!
!
!     3.  Method
!
!     The  LU  decomposition  (or  only  the upper triangular
!     matrix U if the  vector  c  =  L**(-1)*P*b  is  already
!     computed)  is  used  to  obtain an approximation to the
!     solution vector  x  of  the  system  Ax  =  b.  The  LU
!     decomposition  (or only matrix U ) must be available on
!     entry. This means that if a system with a new matrix is
!     solved  then  the  call  of  subroutine  Y12MD  must be
!     preceded by a call of Y12MC and the  contents  of  some
!     parameters and arrays should not be altered between the
!     calls of Y12MC and Y12MD. If a  system  with  the  same
!     matrix  has already been solved, then Y12MD only should
!     be called (but the  contents  of  some  parameters  and
!     arrays  used  in the preceding call of Y12MD should not
!     be altered).
!
!     4.  Parameters of the subroutine
!
!     N     - INTEGER.  On entry N must contain the number of
!             equations in the  system   Ax=b.  Unchanged  on
!             exit.
!
!     A     - REAL  (in  the single precision version Y12MDE)
!             or DOUBLE PRECISION (in  the  double  precision
!             version Y12MDF) array of length NN (see below).
!             On entry array  A  must  contain  the  non-zero
!             elements  of  the  upper  triangular  matrix  U
!             (without the diagonal elements)  and  sometimes
!             also   the   non-zero  elements  of  the  lower
!             triangular matrix L under the  diagonal  (these
!             elements   are  calculated  by  the  subroutine
!             Y12MC, the non-zero elements of L are stored by
!             Y12MC  only when this subroutine is called with
!             IFLAG(5) = 2). The content of the  array  A  is
!             not modified by subroutine Y12MD.
!
!     NN   -  INTEGER. On entry NN must contain the length of
!             the arrays A and SNR.   Restriction:  NN  >  Z.
!             Recommended  value:  2*Z  <  NN  <  3*Z  if the
!             non-zero elements of matrix L will  be  removed
!             during  the  decomposition  (i.e. if subroutine
!             Y12MC is called with IFLAG(5) = 1) and 3*Z < NN
!             < 5*Z if the non-zero elements of matrix L will
!             be stored during the decomposition (i.e. if the
!             subroutine  Y12MC is called with IFLAG(5) = 2).
!             Unchanged on exit.
!
!     B -     REAL (in the single precision  version  Y12MDE)
!             or  DOUBLE  PRECISION  (in the double precision
!             version  Y12MDF)  array   of   length   N.   If
!             subroutine   Y12MC   has  been  called  in  the
!             solution  of  system   Ax   =   b   (i.e.    if
!             IFLAG(5)<3), then on entry array B must contain
!             vector c = L**(-1)*P*b (computed by Y12MC).  If
!             subroutine  Y12MC  has  not  been called in the
!             solution of system Ax = b (i.e  a  system  with
!             the  same matrix has been solved before solving
!             Ax = b and the old LU decomposition is used  by
!             setting  IFLAG(5)  =  3, see below) the array B
!             must contain the right hand side  vector  b  on
!             entry.   The content of the array B is modified
!             by the subroutine Y12MD. On successful exit  an
!             approximation  to  the  true solution x will be
!             found in B.
!
!     PIVOT - REAL (in the single precision  version  Y12MDE)
!             or  DOUBLE  PRECISION  (in the double precision
!             version Y12MDF) array of length N. On entry the
!             diagonal   elements  of  the  upper  triangular
!             matrix U must be stored in array  PIVOT  (these
!             elements  are calculated and stored in PIVOT by
!             subroutine Y12MC).  Unchanged on exit.
!
!     SNR   - INTEGER array of  length  NN  (see  below).  On
!             entry   the  column  numbers  of  the  non-zero
!             elements  of  the  upper  triangular  matrix  U
!             (without  the  column  numbers  of the non-zero
!             elements on the diagonal)  and  sometimes  also
!             the  column numbers of the non-zero elements of
!             the lower  triangular  matrix  L  (without  the
!             column  numbers of the non-zero elements on the
!             diagonal) must be stored in  array  SNR.   This
!             structure  is prepared by Y12MC. The content of
!             array SNR is not  modified  by  the  subroutine
!             Y12MD.
!
!     HA    - INTEGER  two-dimesional  array.   The length of
!             the first dimension of HA is IHA  (see  below).
!             The length of the second dimension of HA is 11.
!             The content of columns 4, 5, 6, 9, 10 and 11 is
!             ignored  by  subroutine  Y12MD.  The content of
!             columns  1,  2,  3,  7  and  8  is  stored   by
!             subroutine  Y12MC  and  should  not  be altered
!             between calls of Y12MC and Y12MD. Unchanged  on
!             exit.
!
!     IHA   - INTEGER.  On  entry IHA must contain the length
!             of   the   first   dimension   of   array   HA.
!             Restriction: IHA .ge. N. Unchanged on exit.
!
!     IFLAG - INTEGER array of length 10. The contents of all
!             locations   in   this  array  except  IFLAG(1),
!             IFLAG(3), IFLAG(4) and IFLAG(5) are ignored  by
!             subroutine  Y12MD. The user should set IFLAG(1)
!             .ge. 0 before the call of  package  Y12M  (i.e.
!             before  the  first call of a subroutine of this
!             package). IFLAG(1) is used in the error  checks
!             by the subroutines and should not be altered by
!             the user between any two  successive  calls  of
!             subroutines  of  the package. IFLAG(1) is equal
!             to  -2  both on entry and  on  successful  exit
!             from  subroutine  Y12MD.  Thus,  the content of
!             IFLAG(1) is not modified by  subroutine  Y12MD.
!             The  same  values  of  IFLAG(3) and IFLAG(4) as
!             those in the preceding call of Y12MC should  be
!             used.  If  subroutine  Y12MC has been called in
!             the solution of the system  Ax  =  b,  then  on
!             entry  IFLAG(5)  must have the same value as on
!             entry of Y12MC. If  subroutine  Y12MC  has  not
!             been  called in the solution of the system Ax =
!             b (i.e. a system with the same matrix has  been
!             solved  before  solving  Ax  =  b  and  the  LU
!             decomposition of matrix A  is  thus  available)
!             then  on entry IFLAG(5) must contain the number
!             3. The content of array IFLAG is  unchanged  on
!             exit.
!
!     IFAIL - Error  diagnostics  parameter.  The  content of
!             parameter  IFAIL  is  modified  by   subroutine
!             Y12MA.  On exit IFAIL = 0 if the subroutine has
!             not detected  any  error.  Positive  values  of
!             IFAIL  on  exit  show  that some error has been
!             detected by the subroutine. Many of  the  error
!             diagnostics  are  common for all subroutines in
!             the package.  Therefore the  error  diagnostics
!             are listed in a separate section, Section 7, of
!             this book. We advise  the  user  to  check  the
!             value of this parameter on exit.
!
!     5.  Error diagnostics
!
!     Error diagnostics are given by positive values  of  the
!     parameter  IFAIL  (see  above).   We advise the user to
!     check carefully the value of this  parameter  on  exit.
!     The error messages are listed in Section 7.
!
!     6.  Auxiliary subroutines
!
!     None.
!
!     7.  Timing
!
!     The  time  taken  depends  on  the  order of the matrix
!     (parameter N) and the number of the  non-zero  elements
!     in the LU-decomposition of matrix A.
!
!     8.  Storage
!
!     There are no internally declared arrays.
!
!     9.  Accuracy
!
!     It  is  difficult  to  evaluate  the  accuracy  of  the
!     computed solution. Large values of parameter  AFLAG(5),
!     the   growth  factor,  indicate  unstable  computations
!     during  the  decomposition  stage.  Small   values   of
!     parameter AFLAG(8), the minimal pivotal element, can be
!     considered as a signal for  numerical  singularity.  We
!     must  emphasize  here  that normally much more reliable
!     evaluations of the accuracy achieved can  be  found  by
!     the  use  of  iterative  refinement, i.e. by the use of
!     subroutine Y12MF. By the use of the  latter  subroutine
!     the  computing time will often be reduced too. However,
!     the storage requirements may sometimes be increased.
!
!     10.  Some remarks
!
!     Remark 1 The content of N, A, SNR, NN, columns 1, 2, 3,
!              7 and 8 of HA,  IHA,  IFLAG(3),  IFLAG(4)  and
!              IFAIL  should  not be altered between calls of
!              Y12MC and Y12MD.
!
!     Remark 2 If   the  LU  decomposition  of  matrix  A  is
!              available (i.e. a system  with  matrix  A  has
!              already  been  solved)  then  Y12MB  and Y12MC
!              should not be called.  The  user  should  only
!              assign the new right hand side vector to array
!              B, set IFLAG(5) = 3 and call Y12MD.
!                =================================
!                Documentation of subroutine Y12MF
!                =================================
!
!     1.  Purpose
!
!     Large  and sparse systems of linear algebraic equations
!     with  real  coefficients  are  solved  by  the  use  of
!     Gaussian   elimination  and  sparse  matrix  technique.
!     Iterative refinement is applied in order to improve the
!     accuracy.
!
!     2.  Calling sequence and declaration of the parameters
!
!     The  subroutine  is  written  in  FORTRAN  and has been
!     extensively tested with the FOR and  FTN  compilers  on
!     the  UNIVAC  1100/82 computer at the Regional Computing
!     Centre at the University of  Copenhagen  (RECKU).  Many
!     examples  have  been run on an IBM 3033 computer at the
!     Northern Europe University Computing Centre (NEUCC) and
!     on  a  CDC Cyber 173 computer at the Regional Computing
!     Centre at the University of  Aarhus  (RECAU).   Only  a
!     single precision version, Y12MFE, of this subroutine is
!     available.  The  call  of  the   subroutine   and   the
!     declarations of the parameters are as follows.
!
!           SUBROUTINE Y12MFE(N, A, SNR, NN, RNR, NN1, A1,
!     SN, NZ,
!          1                  HA, IHA, B, B1, X, Y, AFLAG,
!     IFLAG, IFAIL)
!           REAL A(NN), B(N), B1(N), X(N), Y(N), A1(NZ),
!     AFLAG(11)
!           INTEGER SNR(NN), RNR(NN1), HA(IHA,13), SN(NZ),
!     IFLAG(12)
!           INTEGER N, NN, NN1, NZ, IHA, IFAIL
!
!     This  subroutine  can  be  used on many other computers
!     also. However, some alterations may  be  needed  and/or
!     may ensure greater efficiency of the performance of the
!     subroutine.   For  example,  it  will  be   much   more
!     efficient  to  declare  arrays  SNR,  RNR,  SN  and (if
!     possible)  HA  as  INTEGER*2   arrays   on   some   IBM
!     instalations.  Note too that the subroutine accumulates
!     some inner products in double precision and then rounds
!     them  to  single  precision.  Therefore  if  the use of
!     double precision is very expensive (in comparison  with
!     the  use  of  single  precision) on the computer at the
!     user's disposal then it will not be  efficient  to  use
!     subroutine  Y12MF. The single precision versions of the
!     other subroutines should be used in this situation.
!
!     3.  Method
!
!     Consider the system Ax = b, where matrix A  is  sparse.
!     The  non-zero  elements of matrix A are ordered by rows
!     (subroutine Y12MB is called to perform this  operation)
!     and  then  factorized  into two triangular matrices, an
!     upper triangular matrix U and a unit  lower  triangular
!     matrix  L  (subroutine Y12MC is called to calculate the
!     non-zero elements of U amd L). The factorized system is
!     solved  (subroutine  Y12MD  is  called  to calculate an
!     approximation x1 to  the  solution  vector  x)  and  an
!     attempt  to  improve  the  first  solution by iterative
!     refinement is carried out in the following way:
!
!                      r(k-1) = b - A*x(k-1)
!
!                d(k-1) = Q*U**(-1)*L(-1)*P*r(k-1)
!
!                    x(k)    = X(k-1) + d(k-1)
!
!                     where k = 2, 3, ..., p.
!
!     Normally, this process will improve the accuracy of the
!     first solution x1. If the process is not convergent  or
!     the convegence is very slow, then information about the
!     trouble can be  obtained  by  checking  the  parameters
!     AFLAG(10)  (the  max-norm  of  the last residual vector
!     r(p-1) computed by Y12MF) and AFLAG(9) (the max-norm of
!     the last correction vector, d(p-1), computed by Y12MF).
!     We advise the user to check carefully these  parameters
!     on exit. Large values of these parameters show that the
!     computed solution is not very accurate.
!
!     4.  Parameters of the subroutine
!
!     N     - INTEGER.  On entry N must contain the number of
!             equations in the  system   Ax=b.  Unchanged  on
!             exit.
!
!     A     - REAL  array  of  length  NN (see below). If the
!             coefficient matrix must be factorized (in  this
!             case  IFLAG(5) = 2), then on entry the first NZ
!             locations of array A must contain the  non-zero
!             elements of matrix A, the order of the elements
!             can  be  arbitrary.  If  the  factorization  of
!             matrix  A  is available (i.e. a system with the
!             same matrix has been solved in a previous  call
!             of  Y12MF)  then  on entry array A must contain
!             the non-zero elements of U and L (IFLAG(5) =  3
!             should  also  be  assigned  in  this case). The
!             content of array A is modified  by  Y12MF  when
!             the  coefficient  matrix  has to be factorized.
!             Unchanged on exit when an old LU  factorization
!             is  used.  The content of array A should not be
!             altered between successive calls of  Y12MF  for
!             the solution of systems with the same matrices.
!
!     SNR   - INTEGER  array of the length NN (see below). If
!             the coefficient matrix  A  must  be  factorized
!             then on entry the first NZ positions of array A
!             must contain the column numbers of the non-zero
!             elements  of  matrix A (ordered as the non-zero
!             elements in array A). The content of array  SNR
!             is  modified  by  the  subroutine Y12MF in this
!             case. If the LU factorization of  matrix  A  is
!             available  then on entry array SNR must contain
!             the column numbers of the non-zero elements  of
!             U  and L. The content of array SNR is unchanged
!             on exit in this case. The content of array  SNR
!             should  not be altered between successive calls
!             of Y12MF for the solution of systems  with  the
!             same matrices.
!
!     NN    - INTEGER. On entry NN must contain the length of
!             array A and  SNR.  Restriction:  NN  .ge.  2*Z.
!             Recommended  value:  2*Z   .le.  NN  .le.   3*Z
!             Unchanged on exit.
!
!     RNR   - INTEGER array of length NN1 (see below). If the
!             coefficient  matrix  A must be factorized, then
!             on entry the first NZ positions  of  array  RNR
!             must  contain  the  row numbers of the non-zero
!             elements of matrix A (ordered as  the  non-zero
!             elements  in array A). The content of array RNR
!             is modified by the subroutine in this case.  If
!             the  LU  factorization of matrix A is available
!             then the content of array RNR is ignored by the
!             subroutine Y12MF.
!
!     NN1   - INTEGER.  On  entry NN1 must contain the length
!             of the array RNR.   Restriction:  NN1  .ge.  Z.
!             Recommended  value:  1.5*Z  .le.  NN1 .le. 2*Z.
!             Unchanged on exit.
!
!     A1    - REAL array  of  length  NZ  (see  below).   The
!             content  of  array A1 is modified by subroutine
!             Y12MF  when  the  LU   factorization   of   the
!             coefficient   matrix   A   is   not  available.
!             Subroutine Y12MF copies the first NZ  locations
!             of  array  A  into array A1 in this case (after
!             the internal call of  Y12MB;  this  means  that
!             array  A1 contains the non-zero elements of the
!             coefficient matrix A ordered by rows on  exit).
!             If  the  LU  factorization  of  the coefficient
!             matrix A is available then the content of array
!             A1  is unchanged on exit.  The content of array
!             A1 should  not  be  altered  between  successve
!             calls  of Y12MF in the solution of systems with
!             the same matrices.
!
!     SN    - INTEGER array of length  NZ  (see  below).  The
!             content  of  array SN is modified by subroutine
!             Y12MF  when  the  LU   factorization   of   the
!             coefficient   matrix   A   is   not  available.
!             Subroutine Y12MF copies the first NZ  locations
!             of  array SNR into array SN in this case (after
!             the internal call of  Y12MB;  this  means  that
!             array  SN  contains  the  column numbers of the
!             non-zero elements ordered by rows on exit).  If
!             the  LU  factorization  is  available  then the
!             content of array SN is unchanged on  exit.  The
!             content  of  array  SN  should  not  be altered
!             between  successive  calls  of  Y12MF  in   the
!             solution of systems with the same matrices.
!
!     NZ    - INTEGER. On entry NZ must contain the number of
!             non-zero elements in the coefficient  matrix  A
!             of the system Ax = b. Unchanged on exit.
!
!     HA    - INTEGER  two-dimensional  array.  The length of
!             the first dimension of HA is IHA  (see  below).
!             The length of the second dimension of HA is 13.
!             If a new decomposition  should  be  calculated,
!             then  subroutine Y12MF modifies the contents of
!             HA. The contents on exit of some of the columns
!             of  array HA are described in the documentation
!             of  subroutine  Y12MC.  The  last  two  columns
!             (12'th  and  13'th) contain on exit information
!             about the row starts and the row  ends  in  the
!             original matrix (after the non-zero elements of
!             this matrix have been ordered  by  rows);  this
!             means  that  the  non-zero elements in row i of
!             the coefficient matrix A  are  located  between
!             positions  HA(i,12)  and  HA(i,13)  in array A1
!             (the column numbers of the non-zero elements in
!             row   i  are  also  located  between  positions
!             HA(i,12) and HA(i,13)  in  array  SN).  If  the
!             matrix  of the system has the same structure as
!             the matrix of a system which is already solved,
!             then  subroutine Y12MF will use the information
!             stored in columns 7, 8, 9 and 10 of  HA  during
!             the  previous  call (therefore this information
!             should not be altered). If the LU decomposition
!             of  the  coefficient  matrix is available, then
!             the contents of array HA are unchanged on exit.
!             The  contents  of columns 1, 2, 3, 7, 8, 12 and
!             13 of array HA should not  be  altered  between
!             successive  calls  of  Y12MF in the solution of
!             systems with the same matrix.
!
!     IHA   - INTEGER. On entry IHA must contain  the  length
!             of   the   first   dimension   of   array   HA.
!             Restriction: IHA .ge. N. Unchanged on exit.
!
!     B     - REAL array of length N. On entry the right-hand
!             side vector b of  the  system   Ax=b   must  be
!             stored  in  array  B. The content of array B is
!             modified by subroutine  Y12MF.   On  successful
!             exit  the  components  of  the  last correction
!             vector d(p-1) are stored in array B.
!
!     B1 -    REAL array of length N. The content of array B1
!             is modified by subroutine Y12MF. The right hand
!             side vector of the system Ax = b is  stored  in
!             array B1 on successful exit.
!
!     X -     REAL  array of length N. The content of array X
!             is  modified  by  the  subroutine   Y12MF.   On
!             successful  exit  the corrected solution vector
!             is stored in array X.
!
!     Y -     REAL array of length N. The content of array  Y
!             is  modified by subroutine Y12MF. On successful
!             exit array Y will contain the pivotal  elements
!             (the diagonal elements of matrix U). This means
!             that a small element  (or  small  elements)  in
!             array   Y   on   exit  may  indicate  numerical
!             singularity of the coefficient matrix  A.  Note
!             that  the smallest in absolute value element in
!             array Y is also stored in AFLAG(8), see  below.
!
!     AFLAG - REAL array of length 11.  The  content  of  the
!             array can be described as follows:
!
!             AFLAG(1) -  Stability factor. An element can be
!                         chosen  as  pivotal element only if
!                         this element is larger (in absolute
!                         value)  than  the absolute value of
!                         the  largest  element  in  its  row
!                         divided   by   AFLAG(1).  On  entry
!                         AFLAG(1)  should  contain  a   real
!                         number  larger than 1.0. If this is
!                         not so  then  the  subroutine  sets
!                         AFLAG(1)   =   1.0005.  Recommended
!                         values of AFLAG(1) ranging from 4.0
!                         to  16.0.   Unchanged on exit (when
!                         correctly initialized).
!
!             AFLAG(2) -  Drop-tolerance. An element which in
!                         the  process  of  the  computations
!                         becomes smaller(in absolute  value)
!                         than  the drop-tolerance is removed
!                         from  array  A  (and  its  row  and
!                         column  numbers  are  removed  from
!                         arrays RNR  and  SNR).  Recommended
!                         value:   on entry  AFLAG(2)  should
!                         be in the interval
!                         (a*1.0E-4,a*1.0E-1)  where a is the
!                         magnitude of the elements.  If  the
!                         magnitude    a    of  the  non-zero
!                         elements is  not  known,  then  the
!                         user  should intialize the required
!                         drop-tolerance multiplied  by   -1.
!                         Let a(i) be the maximal element (in
!                         absolute value) in row i (  i=1(1)N
!                         ).   Let  a  be  the minimal number
!                         among all a(i). If a negative value
!                         of   AFLAG(2)   is   assigned,  the
!                         subroutine computes  a  and uses  a
!                         drop-tolerance       equal       to
!                         -AFLAG(2)*a.     In    the    above
!                         considerations  it  is assumed that
!                         the coefficient matrix A is not too
!                         badly  scaled.  If  this  is not so
!                         then perform  row  scaling  of  the
!                         coefficient   matrix  so  that  the
!                         largest elements in the rows are of
!                         magnitude  1 and choose AFLAG(2) in
!                         the    interval    (1.0E-5,1.0E-3).
!                         Unchanged on exit.
!
!             AFLAG(3) -  The   subroutine   will   stop  the
!                         computation when the growth  factor
!                         (parameter   AFLAG(5),  see  below)
!                         becomes larger  than  AFLAG(3).  On
!                         entry  AFLAG(3)  should  contain  a
!                         large    positive    number.     If
!                         AFLAG(3)<1.0E5  then the subroutine
!                         sets  AFLAG(3)=1.0E5.   Recommended
!                         value   AFLAG(3)=1.0E16.  Unchanged
!                         on     exit     (when     correctly
!                         initialized).
!
!             AFLAG(4) -  The   subroutine   will   stop  the
!                         computation when the absolute value
!                         of  a  current  pivotal  element is
!                         smaller   than    AFLAG(4)*AFLAG(6)
!                         (parameter  AFLAG(6)  is  described
!                         below).  On  entry  AFLAG(4)   must
!                         contain    a   small   non-negative
!                         number.  If  AFLAG(4)<0   then  the
!                         subroutine sets
!                         AFLAG(4)=-AFLAG(4).     Recommended
!                         value   AFLAG(4)=1.0E-12. Unchanged
!                         on     exit     (when     correctly
!                         intialized).
!
!             AFLAG(5) -  Growth   factor.   The  content  of
!                         parameter AFLAG(5) is  modified  by
!                         subroutine  Y12MF. After each stage
!                         of   the    Gaussian    elimination
!                         subroutine  Y12MF  sets  AFLAG(5) =
!                         AFLAG(7)/AFLAG(6)       (parameters
!                         AFLAG(6) and AFLAG(7) are described
!                         below). On  exit  large  values  of
!                         parameters  AFLAG(5)  indicate that
!                         an   appreciable   error   in   the
!                         computed  solution  is possible. In
!                         an extreme case, where
!                         AFLAG(5)>AFLAG(3),  the  subroutine
!                         will terminate the computations  in
!                         an attempt to prevent overflow.
!
!             AFLAG(6) -  On  exit  AFLAG(6)  is equal to the
!                         largest element in the  coefficient
!                         matrix  A of the system  Ax=b  (set
!                         by subroutine Y12MB). Unchanged  on
!                         exit.
!
!             AFLAG(7) -  On  exit  the  largest (in absolute
!                         value)  element  found  during  any
!                         stage  of  the  elimination will be
!                         stored in AFLAG(7). The content  of
!                         parameter  AFLAG(7)  is modified by
!                         subroutine Y12MF.
!
!             AFLAG(8) -  On exit the  minimal  (in  absolute
!                         value)   pivotal  element  will  be
!                         stored in AFLAG(8). Small values of
!                         AFLAG(8)     indicate     numerical
!                         singularity  of   the   coefficient
!                         matrix  A.  We  advise  the user to
!                         check this parameter on  exit  from
!                         the calculation very carefully. The
!                         content of  parameter  AFLAG(8)  is
!                         modified by the subroutine Y12MF.
!
!             AFLAG(9) -  The content of AFLAG(9) is modified
!                         by  the   subroutine   Y12MF.   The
!                         max-norm  of  the  last  correction
!                         vector d(p-1)  will  be  stored  in
!                         AFLAG(9)  on  successful  exit.  We
!                         advise the user to check  carefully
!                         this parameter.
!
!             AFLAG(10) - The   content   of   AFLAG(10)   is
!                         modified by the  subroutine  Y12MF.
!                         The  max-norm  of the last residual
!                         vector r(p-1)  will  be  stored  in
!                         AFLAG(10)  on  successful  exit. We
!                         advise the user to check  carefully
!                         this parameter.
!
!             AFLAG(11) - The   content   of   AFLAG(11)   is
!                         modified by the  subroutine  Y12MF.
!                         On  exit AFLAG(11) will contain the
!                         max-norm of the corrected  solution
!                         vector.  If  the value of AFLAG(11)
!                         is  not  to  close  to  zero   then
!                         AFLAG(9)/AFLAG(11)   will  normally
!                         give an excellent evaluation of the
!                         relative   error  in  the  solution
!                         vector.
!
!     IFLAG - INTEGER array of length 12. The content of this
!             array can be described as follows:
!
!             IFLAG(1) - This  parameter  is  used  as a work
!                        space by subroutine Y12MF.
!
!             IFLAG(2) - On entry IFLAG(2) must contain  some
!                        positive  integer smaller than N. We
!                        recommend  IFLAG(2)   .le.   3.   If
!                        IFLAG(3)  = 0 then this parameter is
!                        ignored  by  subroutine  Y12MF.   If
!                        IFLAG((2)  .ge.  0  then the pivotal
!                        search   at   any   stage   of   the
!                        elimination (except possibly some of
!                        the last stages) is carried  out  in
!                        the  IFLAG(2)  rows which have least
!                        number   of    non-zero    elements.
!                        Unchanged on exit.
!
!             IFLAG(3) - On  entry IFLAG(3) must contain 0, 1
!                        or 2.  For  general  pivotal  search
!                        IFLAG(3)  should  be set equal to 1.
!                        If IFLAG(3) = 2 then  only  diagonal
!                        elements of the coefficient matrix A
!                        can be selected as pivotal elements.
!                        If  IFLAG(3)  =  0 then the Gaussian
!                        elimination  will  be  carried   out
!                        without any pivoting.  IFLAG(3)=0 or
!                        IFLAG(3)=2 (i.e. one of the  special
!                        pivotal strategies is to be applied)
!                        should  be   used   very   carefully
!                        because    the   error   diagnostics
!                        algorithm may not trace  all  errors
!                        in this case.  Unchanged on exit.
!
!             IFLAG(4) - On  entry IFLAG(4) must contain 0, 1
!                        or 2.  IFLAG(4)  =  0  is  the  best
!                        choice  when  (i) only one system is
!                        to be solved, (ii) the first  system
!                        of  a  sequence  of systems with the
!                        same  matrix  (Ax   =   b1,   Ax   =
!                        b2,......Ax  =  bp) is to be solved,
!                        (iii) when any system in a  sequence
!                        of  systems  whose  matrices  are of
!                        different structure is to be solved.
!                        IFLAG(4) = 1 is the best choice when
!                        the first system of  a  sequence  of
!                        systems  whose  matrices  are of the
!                        same structure is to be  solved.  In
!                        this  case  IFLAG(4) = 2 can be used
!                        in the solution of any system  after
!                        the first one.  Unchanged on exit.
!
!             IFLAG(5) - If   the  LU  factorization  of  the
!                        coefficient matrix is not available,
!                        then  IFLAG(5)  must  be set to 2 on
!                        entry. If the  LU  factorization  of
!                        the coefficient matrix is available,
!                        then IFLAG(5) must be set  to  3  on
!                        entry.  Unchanged on exit.
!
!             IFLAG(6) - On  successful exit IFLAG(6) will be
!                        equal to  the  number  of  "garbage"
!                        collections in the row ordered list.
!                        If IFLAG(6)  is  large  then  it  is
!                        better  to  choose a larger value of
!                        NN with  next  calls  of  subroutine
!                        Y12MF   with  the  same  or  similar
!                        matrix  A.  This  will  lead  to   a
!                        reduction  in  the  computing  time.
!                        The content of IFLAG(6) is  modified
!                        by the subroutine Y12MF.
!
!             IFLAG(7) - On  successful exit IFLAG(7) will be
!                        equal to  the  number  of  "garbage"
!                        collections  in  the  column ordered
!                        list.  If IFLAG(7) is large then  it
!                        is  better  to choose a larger value
!                        of  NN1  in  the   next   calls   of
!                        subroutine  Y12MF  with  the same or
!                        similar matrix A. This will lead  to
!                        a  reduction  of the computing time.
!                        The content of IFLAG(7) is  modified
!                        by subroutine Y12MF.
!
!             IFLAG(8) - On  successful exit IFLAG(8) will be
!                        equal  to  the  maximal  number   of
!                        non-zero elements kept in array A at
!                        any   stage    of    the    Gaussian
!                        elimination.  If  IFLAG(8)  is  much
!                        smaller than NN (or  NN1)  then  the
!                        length  NN  (or  NN1)  can be chosen
!                        smaller in next calls of  subroutine
!                        Y12MF   with  the  same  or  similar
!                        matrix  A.  This  will  lead  to   a
!                        reduction of the storage needed. The
!                        content of IFLAG(8) is  modified  by
!                        subroutine Y12MF.
!
!             IFLAG(9) - The  minimal  length  NN1  such that
!                        Y12MF  can   solve   systems   whose
!                        matrices  are  of the same structure
!                        without "garbage" collections in the
!                        column ordered list and "movings" of
!                        columns at the  end  of  the  column
!                        ordered  list  is stored in IFLAG(9)
!                        after  the  solution  of  the  first
!                        system   in   the   sequence   (with
!                        IFLAG(4)=1).    The    content    of
!                        IFLAG(9)  is  modified by subroutine
!                        Y12MF when IFLAG(4) = 1 and  ignored
!                        otherwise.
!
!             IFLAG(10)  The  minimal  length  NN  such  that
!                        subroutine Y12MF can  solve  systems
!                        whose   matrices  are  of  the  same
!                        structure     without      "garbage"
!                        collections  in the row ordered list
!                        and "movings" of rows to the end  of
!                        the  row  ordered  list is stored in
!                        IFLAG(10) after the solution of  the
!                        first  system  in the sequence (with
!                        IFLAG(4)  =  1).  The   content   of
!                        IFLAG(10)   is   modified   by   the
!                        subroutine Y12MF when IFLAG(4)  =  1
!                        and ignored otherwise.
!
!             IFLAG(11)  The   maximum   allowed   number  of
!                        iterations  must  be  contained   in
!                        IFLAG(11)  on  entry.   Restriction:
!                        IFLAG(11)>1.   Recommended  value  :
!                        IFLAG(11)< 33.  Unchanged on exit.
!
!             IFLAG(12)  The content of IFLAG(12) is modified
!                        by the subroutine  Y12MF.   On  exit
!                        the  number  of  iterations actually
!                        performed   will   be   stored    in
!                        IFLAG(12).
!
!     IFAIL - Error  diagnostic  parameter.  The  content  of
!             parameter   IFAIL  is  modified  by  subroutine
!             Y12MF.  On exit IFAIL = 0 if the subroutine has
!             not  detected  any  error.  Positive  values of
!             IFAIL on exit show that  some  error  has  been
!             detected  by  the subroutine. Many of the error
!             diagnostics are common for all  subroutines  in
!             the  package.   Therefore the error diagnostics
!             are listed in a separate section, Section 7, of
!             this  book.  We  advise  the  user to check the
!             value of this parameter on exit.
!
!     5.  Error diagnostics
!
!     Error  diagnostics  are  given  by  positive  values of
!     parameter IFAIL (see above).  We  advise  the  user  to
!     check  carefully  the  value of this parameter on exit.
!     The error messages are listed in Section 7.
!
!     6.  Auxiliary subroutines
!
!     Y12MF calls three other subroutines: Y12MB,  Y12MC  and
!     Y12MD.
!
!     7.  Timing
!
!     The  time  taken  depends  on  the  order of the matrix
!     (parameter N), the number of the non-zero  elements  in
!     the matrix (parameter Z), the magnitude of the non-zero
!     elements and their distribution in the matrix.
!
!     8.  Storage
!
!     There are no internally declared arrays.
!
!     9.  Accuracy
!
!     Normally the accuracy achieved can  be  estimated  very
!     well  by  means  of  AFLAG(9)  and  AFLAG(10).  Further
!     information  about  the  accuracy  can   sometimes   be
!     obtained  by  inspection of the growth factor AFLAG(5),
!     and the smallest pivotal element  (in  absolute  value)
!     stored in AFLAG(8).
!
!     10.  Some remarks
!
!     Remark 1 Following the recommendations given in Section
!              1.3.6 the user can write a subroutine (similar
!              to the subroutine Y12MA) where the recommended
!              values of the parameters AFLAG(1) to AFLAG(4),
!              IFLAG(2)  to  IFLAG(5)   and   IFLAG(11)   are
!              initialized.  If this is done then only N, NN,
!              NN1, IHA, A, SNR, RNR and B should be assigned
!              before the call of this new subroutine.
!
!     Remark 2 The  information  stored in some of the arrays
!              should not be altered between successive calls
!              of  subroutine  Y12MF  (see  more  details  in
!              Section 1.3.6).
!                        =================
!                        Error diagnostics
!                        =================
!
!     IFAIL is the error diagnostics parameter.  On exit from
!     each  subroutine  IFAIL  =  0  if  no  error  has  been
!     detected. IFAIL > 0 indicates that some error has  been
!     detected  and  the  computations  have  been terminated
!     immediately after  the  detection  of  the  error.  The
!     errors  corresponding  to the different positive values
!     of IFAIL are listed below.
!
!     IFAIL = 1    The    coefficient   matrix   A   is   not
!                  factorized, i.e. the  call  of  subroutine
!                  Y12MD  was not preceded by a call of Y12MC
!                  during the solution of   Ax=b   or  during
!                  the  solution  of  the  first  system in a
!                  sequence ( Ax1 = b1 , Ax2 = b2,.....,Axp =
!                  bp)  of  systems with the same coefficient
!                  matrix. This will work in all  cases  only
!                  if  the  user  sets IFLAG(1) .ge. 0 before
!                  the call of package Y12M (i.e. before  the
!                  first   call   of  a  subroutine  of  this
!                  package).
!
!     IFAIL = 2    The coefficient matrix A is  not  ordered,
!                  i.e.  the call of subroutine Y12MC was not
!                  preceded by a call  of  Y12MB.  This  will
!                  work  in  all  cases only if the user sets
!                  IFLAG(1) .ge. 0 before the call of package
!                  Y12M  (i.e.  before  the  first  call of a
!                  subroutine of this package).
!
!     IFAIL = 3    A pivotal element abs(a(i,i;j)) < AFLAG(4)
!                  *  AFLAG(6) is selected.  When AFLAG(4) is
!                  sufficiently small this is  an  indication
!                  that the coefficient matrix is numerically
!                  singular.
!
!     IFAIL = 4    AFLAG(5), the  growth  factor,  is  larger
!                  than    AFLAG(3).    When    AFLAG(3)   is
!                  sufficiently large this indicates that the
!                  elements  of the coefficient matrix A grow
!                  so quickly during the  factorization  that
!                  the continuation of the computation is not
!                  justified.  The  choice   of   a   smaller
!                  stability   factor,   AFLAG(1),  may  give
!                  better results in this case.
!
!     IFAIL = 5    The length NN of arrays A and SNR  is  not
!                  sufficient.   Larger  values  of  NN  (and
!                  possibly of NN1) should be used.
!
!     IFAIL = 6    The  length  NN1  of  array  RNR  is   not
!                  sufficient.   Larger  values  of  NN1 (and
!                  possibly of NN) should be used.
!
!     IFAIL = 7    A row without  non-zero  elements  in  its
!                  active    part   is   found   during   the
!                  decomposition.  If   the   drop-tolerance,
!                  AFLAG(2),   is  sufficiently  small,  then
!                  IFAIL = 7 indicates  that  the  matrix  is
!                  numerically  singular. If a large value of
!                  the drop-tolerance AFLAG(2) is used and if
!                  IFAIL = 7  on exit, this is not certain. A
!                  run  with  a  smaller  value  of  AFLAG(2)
!                  and/or  a  careful check of the parameters
!                  AFLAG(8) and AFLAG(5)  is  recommended  in
!                  the latter case.
!
!     IFAIL = 8    A  column without non-zero elements in its
!                  active   part   is   found   during    the
!                  decomposition.   If   the  drop-tolerance,
!                  AFLAG(2),  is  sufficiently  small,   then
!                  IFAIL  =  8  indicates  that the matrix is
!                  numerically singular. If a large value  of
!                  the drop-tolerance AFLAG(2) is used and if
!                  IFAIL = 8  on exit, this is not certain. A
!                  run  with  a  smaller  value  of  AFLAG(2)
!                  and/or a careful check of  the  parameters
!                  AFLAG(8)  and  AFLAG(5)  is recommended in
!                  the latter case.
!
!     IFAIL = 9    A pivotal element  is  missing.  This  may
!                  occur  if  AFLAG(2)  >  0 and IFLAG(4) = 2
!                  (i.e. some system after the first one in a
!                  sequence   of   systems   with   the  same
!                  structure is solved using a positive value
!                  for  the drop-tolerance). The value of the
!                  drop-tolerance   AFLAG(2),    should    be
!                  decreased  and  the  coefficient matrix of
!                  the system refactorized.  This  error  may
!                  also occur when one of the special pivotal
!                  strategies (IFLAG(3)=0 or  IFLAG(3)=2)  is
!                  used  and  the  matrix is not suitable for
!                  such a strategy.
!
!     IFAIL = 10   Subroutine Y12MF is called with IFLAG(5) =
!                  1  (i.e.  with  a  request  to  remove the
!                  non-zero elements of the lower  triangular
!                  matrix    L).     IFLAG(5)=2     must   be
!                  initialized instead of IFLAG(5)=1.
!
!     IFAIL = 11   The coefficient matrix A contains at least
!                  two  elements  in the same position (i,j).
!                  The  input   data   should   be   examined
!                  carefully in this case.
!
!     IFAIL = 12   The number of equations in the system Ax=b
!                  is smaller than 2 (i.e.  N<2).  The  value
!                  of N should be checked.
!
!     IFAIL = 13   The  number  of  non-zero  elements of the
!                  coefficient matrix is  non-positive  (i.e.
!                  Z.le.0  ).   The  value of the parameter Z
!                  (renamed NZ in Y12MF) should be checked.
!
!     IFAIL = 14   The number of  non-zero  elements  in  the
!                  coefficient  matrix  is  smaller  than the
!                  number of equations (i.e. Z  <  N  ).   If
!                  there  is no mistake (i.e. if parameter Z,
!                  renamed NZ in Y12MF, is correctly assigned
!                  on  entry)  then the coefficient matrix is
!                  structurally singular in this case.
!
!     IFAIL = 15   The length IHA of the first  dimension  of
!                  array  HA  is  smaller  than  N.  IHA.ge.N
!                  should be assigned.
!
!     IFAIL = 16   The value of  parameter  IFLAG(4)  is  not
!                  assigned  correctly.  IFLAG(4)  should  be
!                  equal to 0, 1 or 2. See  more  details  in
!                  the description of this parameter.
!
!     IFAIL = 17   A  row  without non-zero elements has been
!                  found in the coefficient matrix A  of  the
!                  system  before the Gaussian elimination is
!                  initiated.  Matrix   A   is   structurally
!                  singular.
!
!     IFAIL = 18   A  column  without  non-zero  elements has
!                  been found in the coefficient matrix A  of
!                  the system before the Gaussian elimination
!                  is initiated.  Matrix  A  is  structurally
!                  singular.
!
!     IFAIL = 19   Parameter  IFLAG(2) is smaller than 1. The
!                  value of IFLAG(2)  should  be  a  positive
!                  integer (IFLAG(2) = 3 is recommended).
!
!     IFAIL = 20   Parameter   IFLAG(3)   is  out  of  range.
!                  IFLAG(3) should be equal to 0, 1 or 2.
!
!     IFAIL = 21   Parameter  IFLAG(5)  is  out   of   range.
!                  IFLAG(5) should be equal to 1, 2 or 3 (but
!                  when IFLAG(5) = 3 Y12MB and  Y12MC  should
!                  not  be  called;  see also the message for
!                  IFAIL = 22 below).
!
!     IFAIL = 22   Either  subroutine  Y12MB  or   subroutine
!                  Y12MC is called with IFLAG(5) = 3. Each of
!                  these subroutines should  be  called  with
!                  IFLAG(5) equal to 1 or 2.
!
!     IFAIL = 23   The    number    of   allowed   iterations
!                  (parameter IFLAG(11) when Y12MF  is  used)
!                  is  smaller  than  2.   IFLAG(11)  .ge.  2
!                  should be assigned.
!
!     IFAIL = 24   At least one element whose  column  number
!                  is  either larger than N or smaller than 1
!                  is found.
!
!     IFAIL = 25   At least one element whose row  number  is
!                  either  larger than N or smaller than 1 is
!                  found.
!                           ==========
!                           References
!                           ==========
!        1.  Bjorck, A. -
!                  "Methods for Sparse  Linear  Least-Squares
!                  Problems".
!                  In: "Sparse Matrix Computations"
!                  (J.Bunch and D.Rose, eds.), pp.177-199.
!                  Academic Press, New York, 1976.
!
!        2.  Clasen, R.J. -
!                  "Techniques  for  Automatic  Tolerance  in
!                  Linear Programming",
!                  Comm. ACM 9, pp. 802-803, 1966.
!
!        3.  Cline, A.K.,  Moler,  C.B.,  Stewart,  G.W.  and
!            Wilkinson, J.H. -
!                  "An estimate for the condition number of a
!                  matrix",
!                  SIAM J. Numer. Anal. 16, 368-375, 1979.
!
!        4.  Dongarra,  J.J.,  Bunch,  J.R.,  Moler, C.B. and
!            Stewart, G.W. -
!                  "LINPACK User's Guide",
!                  SIAM, Philadelphia, 1979.
!
!        5.  Duff, I.S. -
!                  "MA28 - a Set of FORTRAN Subroutines
!                  for Sparse Unsymmetric Matrices",
!                  Report    No.   R8730,   A.E.R.E.,Harwell,
!                  England, 1977.
!
!        6.  Duff, I.S. and Reid, J.K. -
!                  "Some  Design  Features of a Sparse Matrix
!                  Code",
!                  ACM Trans. Math. Software  5,  pp.  18-35,
!                  1979.
!
!        7.  Forsythe, G.E., Malcolm, M.A., and Moler, C.B. -
!                  "Computer   Methods    for    Mathematical
!                  Computations",
!                  Prentice-Hall,   Englewood  Cliffs,  N.J.,
!                  1977.
!
!        8.  Forsythe, G.E. and Moler, C.B. -
!                  "Computer  Solution  of  Linear  Algebraic
!                  Equations",
!                  Prentice-Hall,  Englewood  Cliffs,   N.J.,
!                  1967.
!
!        9.  Gustavson, F.G. -
!                  "Some Basic Techniques for Solving Sparse
!                  Systems of Linear Equations".
!                  In:    "Sparse    Matrices    and    Their
!                  Applications",
!                  (D.J.  Rose and R.A. Willoughby, eds.), pp
!                  41-52,
!                  Plenum Press, New York, 1972.
!
!       10.  Gustavson, F.G. -
!                  "Two  Fast Algorithms for Sparse Matrices:
!                  Multiplication and Permuted
!                  Transposition",
!                  ACM Trans. Math. Software, 4, pp. 250-269,
!                  1978.
!
!       11.  Houbak, N. and Thomsen, P.G. -
!                  "SPARKS  -  a   FORTRAN   Subroutine   for
!                  Solution
!                  of  Large  Systems  of  Stiff  ODE's  with
!                  Sparse Jacobians",
!                  Report  79-02,  Institute  for   Numerical
!                  Analysis,
!                  Technical  University  of Denmark, Lyngby,
!                  Denmark, 1979.
!
!       12.  Moler, C.B. -
!                  "Three  Research  Problems  in   Numerical
!                  Linear Algebra".
!                  In:   "Proceedings   of  the  Symposia  in
!                  Applied Mathematics"
!                  (G.H. Golub and J. Oliger, eds.), Vol. 22,
!                  pp. 1-18,
!                  American Mathematical Society,
!                  Providence, Rhode Island, 1978.
!
!       13.  NAG Library -
!                  Fortran Manual, Mark 7, Vol. 3, 4,
!                  Numerical Algorithms Group,
!                  7 Banbury Road,  Oxford  OX2  6NN,  United
!                  Kingdom.
!
!       14.  Reid, J.K. -
!                  "A  Note  on  the  Stability  of  Gaussian
!                  Elimination",
!                  J.  Inst.  Math.  Appl.,  8,  pp. 374-375,
!                  1971.
!
!       15.  Reid, J.K. -
!                  "Fortran  Subroutines  for Handling Sparse
!                  Linear Programming Bases",
!                  Report R8269, A.E.R.E., Harwell,  England,
!                  1976.
!
!       16.  Schaumburg, K. and Wasniewski, J. -
!                  "Use   of   a   Semiexplicit   Runge-Kutta
!                  Integration  Algorithm  in a Spectroscopic
!                  Problem"
!                  Computers  and  Chemistry  2,  pp.  19-25,
!                  1978.
!
!       17.  Schaumburg, K., Wasniewski, J. and Zlatev, Z. -
!                  "Solution   of    Ordinary    Differential
!                  Equations.   Development of a Semiexplicit
!                  Runge-Kutta Algorithm and Application to a
!                  Spectroscopic Problem",
!                  Computers  and  Chemistry,  3,  pp. 57-63,
!                  1979.
!
!       18.  Schaumburg, K., Wasniewski, J. and Zlatev, Z. -
!                  "The Use of Sparse Matrix Technique in the
!                  Numerical Integration of Stiff Systems  of
!                  Linear Ordinary Differential Equations",
!                  Computers  and  Chemistry,  4,  pp.  1-12,
!                  1980.
!
!       19.  Stewart, G.W. -
!                  "Introduction to Matrix Computations",
!                  Academic Press, New York, 1973.
!
!       20.  Tewarson, R.P. -
!                  "Sparse Matrices",
!                  Academic Press, New York, 1973.
!
!       21.  Thomsen, P.G. -
!                  "Numerical  Solution of Large Systems with
!                  Sparse Jacobians".
!                  In: "Working Papers for  the  1979  SIGNUM
!                  Meeting on Numerical Ordinary Differential
!                  Equations" (R.D. Skeel, ed.),
!                  Computer Science Department,
!                  University  of  Illinois   at   Urbana   -
!                  Champaign, Urbana, Illinois, 1979.
!
!       22.  Wasniewski, J., Zlatev, Z. and Schaumburg, K. -
!                  "A Multibanking  Option  of  an  Iterative
!                  Refinement Subroutine".
!                  In:  "Conference Proceedings and Technical
!                  Papers",
!                  Spring   Conference   of   Univac    Users
!                  Assotiation/Europe, Geneva, 1981.
!
!       23.  Wilkinson, J.H. -
!                  "Rounding Errors in Algebraic Processes",
!                  Prentice-Hall,  Englewood  Cliffs,   N.J.,
!                  1963.
!
!       24.  Wilkinson, J.H. -
!                  "The Algebraic Eigenvalue Problem",
!                  Oxford University Press, London, 1965.
!
!       25.  Wilkinson, J.H and Reinsch, C. -
!                  "Handbook for Automatic Computation",
!                  Vol II, Linear Algebra, pp. 50-56,
!                  Springer, Heidelberg, 1971.
!
!       26.  Wolfe, P. -
!                  "Error   in   the   Solution   of   Linear
!                  Programming Problems",
!                  In:   "Error   in   Digital   Computation"
!                  (L.B.Rall, ed.), Vol 2, pp. 271-284,
!                  Wiley, New York, 1965.
!
!       27.  Zlatev, Z. -
!                  "Use   of   Iterative  Refinement  in  the
!                  Solution of Sparse Linear Systems",
!                  Report 1/79, Institute of Mathematics  and
!                  Statistics,   The   Royal  Veterinary  and
!                  Agricultural University,
!                  Copenhagen, Denmark, 1979
!                  (to appear in SIAM J. Numer. Anal.).
!
!       28.  Zlatev, Z. -
!                  "On  Some  Pivotal  Strategies in Gaussian
!                  Elimination by Sparse Technique",
!                  SIAM J. Numer. Anal. 17, pp. 18-30,  1980.
!
!       29.  Zlatev, Z. -
!                  "On Solving Some Large Linear Problems  by
!                  Direct Methods",
!                  Report   111,   Department   of   Computer
!                  Science,  University  of  Aarhus,  Aarhus,
!                  Denmark, 1980.
!
!       30.  Zlatev, Z. -
!                  "Modified  Diagonally Implicit Runge-Kutta
!                  Methods",
!                  Report No.  112,  Department  of  Computer
!                  Science,  University  of  Aarhus,  Aarhus,
!                  Denmark, 1980
!                  (to appear in SIAM Journal  on  Scientific
!                  and Statistical Computing).
!
!       31.  Zlatev, Z. -
!                  "Comparison  of  Two Pivotal Strategies in
!                  Sparse Plane Rotations",
!                  Report   122,   Department   of   Computer
!                  Science,  University  of  Aarhus,  Aarhus,
!                  Denmark, 1980.
!                  (to appear in  Computers  and  Mathematics
!                  with Applications).
!
!       32.  Zlatev, Z., Barker, V.A. and Thomsen, P.G. -
!                  "SSLEST -  a  FORTRAN  IV  Subroutine  for
!                  Solving Sparse Systems of Linear Equations
!                  (USER's GUIDE)",
!                  Report  78-01,  Institute  for   Numerical
!                  Analysis,
!                  Technical  University  of Denmark, Lyngby,
!                  Denmark, 1978.
!
!       33.  Zlatev, Z. and Nielsen, H.B. -
!                  "Least  - Squares Solution of Large Linear
!                  Problems".
!                  In:  "Symposium i Anvendt Statistik 1980"
!                  (A. Hoskuldsson, K.  Conradsen,  B.  Sloth
!                  Jensen and K.Esbensen, eds.), pp. 17-52.
!                  NEUCC, Technical University of Denmark,
!                  Lyngby, Denmark, 1980.
!
!       34.  Zlatev, Z., Schaumburg, K. and Wasniewski, J. -
!                  "Implementation of an Iterative Refinement
!                  Option in a  Code  for  Large  and  Sparse
!                  Systems".
!                  Computers  and  Chemistry,  4,  pp. 87-99,
!                  1980.
!
!       35.  Zlatev, Z. and Thomsen., P.G. -
!                  "ST  -  a  FORTRAN  IV  Subroutine for the
!                  Solution  of  Large  Systems   of   Linear
!                  Algebraic Equations with Real Coefficients
!                  by Use of Sparse Technique",
!                  Report  76-05,  Institute  for   Numerical
!                  Analysis,
!                  Technical  University  of Denmark, Lyngby,
!                  Denmark, 1976.
!
!       36.  Zlatev, Z. and Thomsen, P.G. -
!                  "An   Algorithm   for   the   Solution  of
!                  Parabolic Partial  Differential  Equations
!                  Based on Finite Element Discretization",
!                  Report   77-09,  Institute  for  Numerical
!                  Analysis,
!                  Technical University of  Denmark,  Lyngby,
!                  Denmark, 1977.
!
!       37.  Zlatev, Z. and Thomsen, P.G. -
!                  "Application of  Backward  Differentiation
!                  Methods  to the Finite Element Solution of
!                  Time Dependent Problems",
!                  International   Journal   for    Numerical
!                  Methods  in  Engineering,  14,  pp. 1051 -
!                  1061, 1979.
!
!       38.  Zlatev, Z., Wasniewski, J. and Schaumburg, K. -
!                  "Comparison  of Two Algorithms for Solving
!                  Large Linear Systems".
!                  Report No 80/9,
!                  Regional   Computing   Centre    at    the
!                  University  of Copenhagen, Vermundsgade 5,
!                  DK-2100 Copenhagen, Denmark, 1980.
!
!       39.  Zlatev, Z., Wasniewski, J. and Schaumburg, K. -
!                  "Classification of the Systems of Ordinary
!                  Differential   Equations   and   Practical
!                  Aspects in the  Numerical  Integration  of
!                  Large Systems",
!                  Computers  and  Chemistry,  4,  pp. 13-18,
!                  1980.
!
!       40.  Zlatev, Z., Wasniewski, J., Schaumburg, K. -
!                  "A Testing Scheme For Subroutines  Solving
!                  Large Linear Problems",
!                  Report No 81/1,
!                  Regional    Computing    Centre   at   the
!                  University of Copenhagen, Vermundsgade  5,
!                  DK-2100 Copenhagen, Denmark, 1981
!                  (to  appear in Computers and Chemistry, 5,
!                  1981).
! Acknowledge-To: <NEUJW@vm.uni-c.dk>


    use const
    implicit none

    real(kind=dp) :: real_flag(8)
    integer :: int_flag(10)
    integer :: error_flag
    integer                      :: vec_size1
    integer                      :: vec_size2
    integer                      :: vec_size3
    real(kind = dp), allocatable :: pivot(:)
    integer, allocatable :: ha(:,:)


    integer, parameter :: dim2ha = 11

  contains
    subroutine solve_sys(n, z, nn, nn1, iha, a, snr, rnr,b, pivot,ha,aflag,iflag,ifail) !y12maf
      implicit none
      integer              :: n
      !< größe des Gleichunssystems
      integer              :: z
      !< Anzahl der Nicht-Null Elemente in der Matrix
      integer              :: nn
      !< Länge des Matrix Vektor, vorgeschlagener Wert: 2* z
      integer              :: nn1
      !< Länge des Row-Vektors, kann etwas kleiner als nn aber größer gleich z sein.
      integer              :: iha
      !< Länge des temp Lösungsvektors, vorgeschlagener wert n*2
      real(kind=dp)        :: a(nn)
      integer              :: snr(nn)
      !< zeile
      integer              :: rnr(nn1)
      !< spalte
      real(kind=dp)        :: pivot(n)
      integer              :: ha(iha,11)
      real(kind=dp)        :: aflag(8)
      integer              :: iflag(10)
      real(kind=dp)        :: b(n)
      integer              :: ifail

      aflag(1)=16.0E0_dp
      aflag(2)=1.0E-12_dp
      aflag(3)=1.0E+16_dp
      aflag(4)=1.0E-12_dp
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=1

      call y12mbf(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

      if(ifail.ne.0) then
         write(*,*) "ERROR DURING Y12MBF",ifail
         stop 1
      end if

      call y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,ifail)

      if(ifail.ne.0)then
         write(*,*) "ERROR DURING Y12MCF",ifail
         stop 1
      end if

      call y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)

      end subroutine !y12maf
!
!
!  the non-zero elements of a sparse matrix a are prepared  in order to
!  solve the system ax=b by use of sparse matrix technique/
!
!
      subroutine y12mbf(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag,iflag,ifail)
      implicit none

      integer :: n
      integer :: z
      integer :: nn1
      integer :: nn
      integer :: iha
      real(kind=dp) :: a(nn)
      integer :: snr(nn)
      integer :: rnr(nn1)
      integer :: ha(iha,11)
      real(kind=dp) :: aflag(8)
      integer :: iflag(10)
      integer :: ifail

      integer :: mode,i,j,index,l1,l2,l3,l4,l5,r
      real(kind=dp) :: gt1,t


      mode=iflag(4)
      ifail=0
      if (n.lt.2) ifail=12
      if (z.le.0) ifail=13
      if (nn.lt.2*z) ifail=5
      if (nn1.lt.z)  ifail=6
      if (ifail.eq.0.and.n.gt.z) ifail=14
      if (iha.lt.n)  ifail=15
      if (mode.lt.0) ifail=16
      if (mode.gt.2) ifail=16
      if (ifail.ne.0) goto 22
      gt1=0.0d0
      do i=1,n
         ha(i,2)=0
         ha(i,3)=0
         ha(i,6)=0
      end do
!
!  find the number of the non-zero elements in each row and column;move
!  the non-zero elements in the end of the arrays a and snr;find the
!  largest non-zero element in a(in absolute value).
!
      do i=1,z
         t = abs(a(i))
         l3 = rnr(i)
         l4=snr(i)
         if (l4.gt.n.or.l4.lt.1) ifail=24
         if (l3.gt.n.or.l3.lt.1) ifail=25
         ha(l3,3)=ha(l3,3)+1
         ha(l4,6)=ha(l4,6)+1
         if(t.gt.gt1)gt1=t
         a(z+i)=a(i)
         snr(z+i)=snr(i)
      end do

      if(ifail.gt.0) goto 22
!
!  store the information of the row starts(in ha(i,1))and of the column
!  starts(in ha(i,4)).
!
      l1=1
      l2=1
      do i=1,n
         l3=ha(i,3)
         l4=ha(i,6)
         if (l3.gt.0) goto 21
         ifail=17
         goto 22
21       if(l4.gt.0) goto 23
         ifail=18
         goto 22
23       if(mode.eq.2) goto 30
         ha(i,9)=l3
         ha(i,10)=l4
         ha(i,11)=0
         ha(l3,2)=ha(l3,2)+1
         ha(i,5)=l3
30       ha(i,1)=l1
         ha(i,4)=l2
         l1=l1+l3
         l2=l2+l4
         ha(i,3)=0
         ha(i,6)=0
      end do
!
!  store the non-zero elements of matrix a(ordered in rows) in the
!  first z locations of the array a.do the same for their column numbers
!
      do i=1,z
         l1=z+i
         l3=rnr(i)
         l2=ha(l3,1)+ha(l3,3)
         a(l2)=a(l1)
         snr(l2)=snr(l1)
         ha(l3,3)=ha(l3,3)+1
      end do
!
!  store the row numbers of the non-zero elements ordered by columns in
!  the first z locations of the array rnr. store information about row
!  ends(in ha(i,3)).
!
      l4=1
      do 70 i=1,n
         if(mode.eq.2) goto 60
         if(ha(i,2).eq.0) goto 60
         ha(i,11)=l4
         l4=l4+ha(i,2)
         ha(i,2)=ha(i,11)
60       ha(i,3)=ha(i,1)+ha(i,3)-1
         l1=ha(i,1)
         l2=ha(i,3)
         do 70 j=l1,l2
         l3=snr(j)
         r=ha(l3,6)
         index=ha(l3,4)+r
         rnr(index)=i
         if(r.eq.0) goto 70
         if(j.eq.l1) goto 70
         if(rnr(index-1).ne.i) goto 70
         ifail=11
         goto 22
70       ha(l3,6)=r+1

      do 90 i=1,n
      if(mode.eq.2) goto 80
      l3=ha(i,5)
      l5=ha(l3,2)
      ha(l5,8)=i
      ha(i,7)=l5
      ha(l3,2)=ha(l3,2)+1
80    continue
90    ha(i,6)=ha(i,4)+ha(i,6)-1
      aflag(6)=gt1
      iflag(6)=0
      iflag(7)=0
      iflag(8)=z
      iflag(1)=-1
22    return
      end subroutine y12mbf







      subroutine y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha, aflag,iflag,ifail)
!
!  systens of linear equations are solved by use of sparse matrix tech-
!  nique and by gaussian elimination.
!
      !implicit double precision(a-b,g,p,t-y),integer(c,f,h-n,r-s,z)
      integer :: n,z,nn,nn1,iha,ifail

      real(kind=dp) ::  a(nn),b(n),pivot(n),aflag(8)
!
!  information which is necessary to begin the elimination is stored.
!
      integer :: snr(nn),rnr(nn1),ha(iha,11), iflag(10)

!!!!! local variables
      integer :: c1,c2,c3,cr1,cr2,cr3,cr4,i,k,kk,i1,index,j,jj
      integer :: l,l1,l2,l3,l4,l5,l6,l7,ll,lfc,lfr
      integer :: nr,n7,n8,zz
      integer :: r,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,rcoll
      integer :: rr,rr1,rr2,rr3,rr4,rpivot,rrow,slut
      real(kind=dp) :: grmin,t,td,td1,tol1,tol2,tol3,u,v


      ifail=0
      if(iflag(1).ne.-1)ifail=2
      if(aflag(1).lt.1.0E0_dp) aflag(1)=1.0005E0_dp
      if(aflag(3).lt.1.0E+5_dp) aflag(3)=1.0E+5_dp
      if(aflag(4).lt.0.0E0_dp) aflag(4)=-aflag(4)
      if(iflag(2).lt.1) ifail=19
      if(iflag(3).lt.0.or.iflag(3).gt.2) ifail=20
      if(iflag(5).lt.1.or.iflag(5).gt.3) ifail=21
      if(iflag(5).eq.3) ifail=22
      if(ifail.gt.0) goto 1110
      snr(z+1)=0
      rnr(z+1)=0
      n8=n+1
      n7=n-1
      u=aflag(1)
      grmin=aflag(4)*aflag(6)
!
!  use the information about fill-ins if it is possible.
!
      zz=z
      nr=n*n
      if(iflag(4).ne.2) goto 100
      if(iflag(10).gt.nn) goto 50
      l1=iflag(10)
      l5=l1+1
      if(l5.le.nn)snr(l5)=0
      do 40 i=1,n
      l=n8-i
      l2=ha(l,3)+1
      l3=l2-ha(l,1)
      do 10 j=1,l3
      snr(l5-j)=snr(l2-j)
   10 a(l5-j)=a(l2-j)
      ha(l,3)=l1
      ha(l,1)=l5-l3
      l6=l1-l3
      l5=l5-ha(l,9)
      if(l5.gt.l6) goto 30
      do 20 j=l5,l6
   20 snr(j)=0
   30 continue
   40 l1=l5-1
   50 if(iflag(9).gt.nn1) goto 100
      l2=iflag(9)
      l5=l2+1
      if(l5.le.nn1)rnr(l5)=0
      do 90 i=1,n
      l=n8-i
      l1=ha(l,6)+1
      l4=l1-ha(l,4)
      do 60 j=1,l4
   60 rnr(l5-j)=rnr(l1-j)
      ha(l,4)=l5-l4
      ha(l,6)=l2
      l6=l2-l4
      l5=l5-ha(l,10)
      if(l5.gt.l6) goto 80
      do 70 j=l5,l6
   70 rnr(j)=0
   80 continue
   90 l2=l5-1
  100 r4=ha(n,3)
      r5=ha(n,6)
      aflag(7)=aflag(6)
      aflag(8)=aflag(6)
      do 110 i=1,n
      pivot(i)=0.0E0_dp
      ha(i,2)=ha(i,1)
  110 ha(i,5)=ha(i,4)
      index=ha(n,8)
!
!  start of gaussian elimination.
!
      slut=ha(index,3)-ha(index,2)+1
      do 950 i=1,n7
      rr3=ha(i,2)
      rr4=ha(i,3)
      c1=ha(i,4)
      cr4=ha(i,6)
      if(iflag(3).eq.0) goto 350
      if(iflag(4).ne.2) goto 120
      rrow=ha(i,7)
      rcoll=ha(i,8)
      goto 220
  120 l4=ha(i,8)
      if(iflag(3).eq.1) goto 130
      rrow=l4
      rcoll=rrow
      rpivot=i
      goto 170
  130 r=nr
      v=0.0E0_dp
      index=iflag(2)
      do 160 kk=1,index
      l1=i-1+kk
      if(l1.gt.n) goto 170
      j=ha(l1,8)
      r7=ha(j,2)
      r8=ha(j,3)
      r9=r8-r7
      t=0.0E0_dp
      do 140 k=r7,r8
      td=dabs(a(k))
  140 if(t.lt.td)t=td
      t=t/u
      do 160 k=r7,r8
      td=dabs(a(k))
      if(td.lt.t) goto 150
      r6=snr(k)
      r3=r9*(ha(r6,6)-ha(r6,5))
      if(r3.gt.r) goto 150
      if(r3.lt.r) goto 151
      if(v.ge.td) goto 150
  151 v=td
      rrow=j
      rcoll=r6
      r=r3
      rpivot=l1
  150 continue
  160 continue
  170 r3=ha(rcoll,10)
      ha(rcoll,10)=ha(i,10)
      ha(i,10)=r3
      r3=ha(rrow,9)
      ha(rrow,9)=ha(i,9)
!
!  remove the pivot row of the list where the rows are ordered by
!  increasing numbers of non-zero elements.
!
      ha(i,9)=r3
      l1=0
      l=i
      l2=ha(l4,3)-ha(l4,2)+1
  180 l=l+1
      if(l2.gt.l1) ha(l2,11)=l
      if(l.gt.n) goto 190
      l5=ha(l,8)
      l3=ha(l5,3)-ha(l5,2)+1
      if(rpivot.lt.l) goto 190
      ha(l4,7)=l
      ha(l,8)=l4
      l4=l5
      l1=l2
      l2=l3
      l3=n8
      goto 180
  190 if(l2.eq.l1) goto 200
      if(l3.eq.l2) goto 200
      ha(l2,11)=0
  200 l5=ha(i,7)
      if(rrow.eq.i) goto 210
      ha(l5,8)=rrow
      ha(rrow,7)=l5
  210 ha(i,7)=rrow
!
!  row interchanges.
!
      ha(i,8)=rcoll
  220 if(rrow.eq.i) goto 290
      t=b(rrow)
      b(rrow)=b(i)
      b(i)=t
      do 250 j=rr3,rr4
      l1=snr(j)
      r=ha(l1,5)-1
      r10=ha(l1,6)
  240 r=r+1
      if(rnr(r).ne.i) goto 240
      rnr(r)=rnr(r10)
  250 rnr(r10)=rrow
      rr3=ha(rrow,2)
      rr4=ha(rrow,3)
      do 270 j=rr3,rr4
      l1=snr(j)
      r=ha(l1,5)-1
  260 r=r+1
      if(rnr(r).ne.rrow) goto 260
  270 rnr(r)=i
      do 280 j=1,3
      r3=ha(rrow,j)
      ha(rrow,j)=ha(i,j)
!
!  column interchanges.
!
  280 ha(i,j)=r3
  290 if(rcoll.eq.i) goto 350
      do 310 j=c1,cr4
      l1=rnr(j)
      r=ha(l1,2)-1
      r10=ha(l1,3)
  300 r=r+1
      if(snr(r).ne.i) goto 300
      t=a(r10)
      a(r10)=a(r)
      a(r)=t
      snr(r)=snr(r10)
  310 snr(r10)=rcoll
      c1=ha(rcoll,4)
      cr4=ha(rcoll,6)
      do 330 j=c1,cr4
      l1=rnr(j)
      r=ha(l1,2)-1
  320 r=r+1
      if(snr(r).ne.rcoll) goto 320
  330 snr(r)=i
      do 340 j=4,6
      r3=ha(rcoll,j)
      ha(rcoll,j)=ha(i,j)
!
! end of the interchanges.
! the row ordered list and the column ordered list are prepared to
! begin step i of the elimination.
!
  340 ha(i,j)=r3
  350 r9=rr4-rr3
      do 360 rr=rr3,rr4
      if(snr(rr).eq.i) goto 370
  360 continue
      ifail=9
      goto 1110
  370 v=a(rr)
      pivot(i)=v
      td=dabs(v)
      if(td.lt.aflag(8))aflag(8)=td
      if(td.ge.grmin) goto 380
      ifail=3
      goto 1110
  380 r2=ha(i,1)
      a(rr)=a(rr3)
      snr(rr)=snr(rr3)
      a(rr3)=a(r2)
      snr(rr3)=snr(r2)
      snr(r2)=0
      z=z-1
      rr3=rr3+1
      ha(i,2)=rr3
      ha(i,1)=r2+1
      cr3=ha(i,5)
      if(r9.le.0) goto 431
      do 430 j=rr3,rr4
      index=snr(j)
  430 pivot(index)=a(j)
  431 r7=cr4-cr3+1
      do 880 k=1,r7
      r1=rnr(cr3-1+k)
      if(r1.eq.i) goto 870
      i1=ha(r1,1)
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      l2=rr2-rr1+1
      l=rr1-1
  390 l=l+1
      if(snr(l).ne.i) goto 390
      t=a(l)/v
      if(iflag(5).eq.2) goto 400
      a(l)=a(i1)
      snr(l)=snr(i1)
      snr(i1)=0
      i1=i1+1
      ha(r1,1)=i1
      z=z-1
      goto 410
  400 a(l)=a(rr1)
      a(rr1)=t
      r3=snr(rr1)
      snr(rr1)=snr(l)
      snr(l)=r3
  410 rr1=rr1+1
      ha(r1,2)=rr1
      b(r1)=b(r1)-b(i)*t
      if(r9.le.0) goto 669
      r=rr1
      if(r.gt.rr2) goto 470
      do 460 l=r,rr2
      l1=snr(l)
      td=pivot(l1)
      if(td.eq.0.0E0_dp) goto 450
      pivot(l1) = 0.0E0_dp
      td=a(l)-td*t
      a(l)=td
      td1=dabs(td)
      if(td1.gt.aflag(7)) aflag(7)=td1
!
!  too small element is created.remove it from the lists.
!
      if(td1.gt.aflag(2)) goto 450
      z=z-1
      a(l)=a(rr1)
      snr(l)=snr(rr1)
      a(rr1)=a(i1)
      snr(rr1)=snr(i1)
      snr(i1)=0
      rr1=rr1+1
      i1=i1+1
      ha(r1,2)=rr1
      ha(r1,1)=i1
      r3=ha(l1,5)
      r2=r3-1
      l4=ha(l1,4)
      l5=rnr(l4)
      l6=rnr(r3)
  440 r2=r2+1
      if(rnr(r2).ne.r1) goto 440
      rnr(r2)=l6
      rnr(r3)=l5
      rnr(l4)=0
      ha(l1,5)=r3+1
      ha(l1,4)=l4+1
  450 continue
  460 continue
  470 continue
      do 750 j=1,r9
      r=rr3-1+j
      r2=snr(r)
      tol2=pivot(r2)
      pivot(r2)=a(r)
      if(tol2.eq.0.0E0_dp) goto 740
      tol3=-tol2*t
      tol1=dabs(tol3)
      if(tol1.lt.aflag(2)) goto 740
      c2=ha(r2,4)
      cr2=ha(r2,6)
      cr1=ha(r2,5)
      lfr=rr2-i1+2
      lfc=cr2-c2+2
      if(iflag(4).ne.1) goto 480
      if(lfr.gt.ha(r1,9))ha(r1,9)=lfr
      if(lfc.gt.ha(r2,10))ha(r2,10)=lfc
  480 if(i1.eq.1) goto 490
      if(snr(i1-1).eq.0) goto 600
  490 if(rr2.eq.nn) goto 500
      if(snr(rr2+1).eq.0) goto 580
!
!  collection in row ordered list.
!
  500 r10=nn-lfr
      if(r10.ge.r4) goto 560
      iflag(6)=iflag(6)+1
      do 520 jj=1,n
         l1=ha(jj,3)
         if(l1.lt.ha(jj,1)) goto 510
         ha(jj,3)=snr(l1)
         snr(l1)=-jj
  510 continue
  520 continue
      l3=0
      l4=1
      do 550 jj=1,r4
         if(snr(jj).eq.0) goto 540
         l3=l3+1
         if(snr(jj).gt.0) goto 530
         l5=-snr(jj)
         snr(jj)=ha(l5,3)
         ha(l5,3)=l3
         l6=l4+ha(l5,2)-ha(l5,1)
         ha(l5,2)=l6
         ha(l5,1)=l4
         l4=l3+1
     530 a(l3)=a(jj)
         snr(l3)=snr(jj)
     540 continue
  550 continue
      r4=l3
      snr(l3+1)=0
      rr3=ha(i,2)
      rr4=ha(i,3)
      i1=ha(r1,1)
      rr1=ha(r1,2)
      r=rr3-1+j
      !write(*,*) r10,r4
      if(r10.ge.r4) goto 560
      ifail=5
!
! fill-in takes place in the row ordered list.
!
      goto 1110
  560 r8=lfr-1
      rr2=r4+lfr
      if(r8.le.0) goto 579
      l3=i1-1
      do 570 ll=1,r8
      l4=r4+ll
      l5=l3+ll
      a(l4)=a(l5)
      snr(l4)=snr(l5)
  570 snr(l5)=0
  579 rr1=r4+rr1-i1+1
      ha(r1,3)=rr2
      ha(r1,2)=rr1
      i1=r4+1
      ha(r1,1)=i1
      l1=rr2
      goto 590
  580 rr2=rr2+1
      ha(r1,3)=rr2
      l1=rr2
      if(rr2.le.r4) goto 610
  590 r4=rr2
      if(r4.lt.nn)snr(r4+1)=0
      goto 610
  600 rr1=rr1-1
      i1=i1-1
      ha(r1,1)=i1
      ha(r1,2)=rr1
      l1=rr1
      snr(i1)=snr(l1)
      a(i1)=a(l1)
  610 a(l1)=tol3
      snr(l1)=snr(r)
      td=dabs(a(l1))
      if(td.gt.aflag(7))aflag(7)=td
      z=z+1
      if(iflag(8).lt.z) iflag(8)=z
      if(c2.eq.1) goto 620
      if(rnr(c2-1).eq.0) goto 720
  620 if(cr2.eq.nn1) goto 630
      if(rnr(cr2+1).eq.0) goto 700
!
!  collection in column ordered list.
!
  630 r10=nn1-lfc
      if(r10.ge.r5) goto 680
      iflag(7)=iflag(7)+1
      do 640 jj=i,n
      l1=ha(jj,6)
      ha(jj,6)=rnr(l1)
  640 rnr(l1)=-jj
      l3=0
      l4=1
      do 670 jj=1,r5
      if(rnr(jj).eq.0) goto 660
      l3=l3+1
      if(rnr(jj).gt.0) goto 650
      l5=-rnr(jj)
      rnr(jj)=ha(l5,6)
      ha(l5,6)=l3
      l6=l4+ha(l5,5)-ha(l5,4)
      ha(l5,5)=l6
      ha(l5,4)=l4
      l4=l3+1
  650 rnr(l3)=rnr(jj)
  660 continue
  670 continue
      r5=l3
      rnr(r5+1)=0
      c2=ha(r2,4)
      cr3=ha(i,5)
      cr4=ha(i,6)
      cr1=ha(r2,5)
      if(r10.ge.r5) goto 680
      ifail=6
!
! fill-in takes place in the column ordered list.
!
      goto 1110
  680 r8=lfc-1
      cr2=r5+lfc
      if(r8.le.0) goto 699
      l3=c2-1
      do 690 l=1,r8
      l4=r5+l
      l5=l3+l
      rnr(l4)=rnr(l5)
  690 rnr(l5)=0
  699 cr1=r5+cr1-c2+1
      c2=r5+1
      ha(r2,6)=cr2
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr2
      goto 710
  700 cr2=cr2+1
      ha(r2,6)=cr2
      r=cr2
      if(cr2.le.r5) goto 730
  710 r5=cr2
      if(r5.lt.nn1)rnr(r5+1)=0
      goto 730
  720 cr1=cr1-1
      c2=c2-1
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr1
      rnr(c2)=rnr(r)
  730 rnr(r)=r1
  740 continue
  750 continue
  669 if(rr1.le.rr2) goto 760
      ifail=7
!
!  update the information in the list where the rows are ordered by
!  increasing numbers of the non-zero elements.
!
      goto 1110
  760 if(iflag(4).eq.2) goto 870
      if(iflag(3).eq.0) goto 870
      l1=rr2-rr1+1
      if(l1.eq.l2) goto 870
      l6=ha(r1,7)
      l4=ha(l2,11)
      if(l1.gt.l2) goto 820
      if(l6.gt.l4) goto 780
      if(l4.eq.n) goto 770
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2) goto 790
  770 ha(l2,11)=0
      goto 800
  780 l5=ha(l4,8)
      l3=ha(l6,8)
      ha(l4,8)=l3
      ha(l6,8)=l5
      ha(l5,7)=l6
      ha(l3,7)=l4
      l6=l4
  790 ha(l2,11)=l4+1
  800 if(l4.eq.i+1) goto 810
      l=ha(l6-1,8)
      l2=ha(l,3)-ha(l,2)+1
      l4=ha(l2,11)
      if(l1.lt.l2) goto 780
  810 if(l1.ne.l2)ha(l1,11)=l6
      goto 870
  820 if(l6.gt.l4) goto 840
      if(l4.eq.n) goto 830
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2) goto 840
  830 ha(l2,11)=0
  840 l2=l2+1
      if(l2.le.slut) goto 850
      l3=n
      slut=l1
      l2=l1
      goto 860
  850 l3=ha(l2,11)-1
      if(l3.eq.-1) goto 840
      if(l2.gt.l1)l2=l1
  860 ha(l2,11)=l3
      l4=ha(l3,8)
      l7=ha(l6,8)
      ha(l3,8)=l7
      ha(l6,8)=l4
      ha(l7,7)=l3
      ha(l4,7)=l6
      l6=l3
      if(l2.lt.l1) goto 840
  870 continue
  880 continue
      if(r9.le.0) goto 882
      do 881 j=rr3,rr4
      index=snr(j)
  881 pivot(index)=0.0d0
  882 continue
      cr3=ha(i,4)
      do 890 j=cr3,cr4
  890 rnr(j)=0
      if(r9.le.0) goto 930
      l2=ha(i,2)-1
      do 920 ll=1,r9
      r=snr(l2+ll)
      r1=ha(r,5)
      r2=ha(r,6)
      if(r2.gt.r1) goto 900
      ifail=8
      goto 1110
  900 ha(r,5)=r1+1
      r3=r1-1
  910 r3=r3+1
      if(rnr(r3).ne.i) goto 910
      rnr(r3)=rnr(r1)
  920 rnr(r1)=i
  930 aflag(5)=aflag(7)/aflag(6)
      if(aflag(5).lt.aflag(3)) goto 940
      ifail=4
      goto 1110
  940 continue
!
!  preparation to begin the back substitution.
!
  950 continue
      index=ha(n,2)
      pivot(n)=a(index)
      a(index)=0.0d0
      td=dabs(pivot(n))
      if(td.gt.aflag(7))aflag(7)=td
      if(td.lt.aflag(8))aflag(8)=td
      if(td.gt.grmin) goto 960
      ifail=3
      goto 1110
  960 if(iflag(4).ne.1) goto 1060
      iflag(10)=ha(n,9)
      iflag(9)=ha(n,10)
      do 990 i=1,n7
      r1=n-i
      iflag(10)=iflag(10)+ha(r1,9)
      iflag(9)=iflag(9)+ha(r1,10)
      if(iflag(3).eq.0) goto 980
      do 970 j=9,10
      r2=ha(r1,j-2)
      r6=ha(r2,j)
      ha(r2,j)=ha(r1,j)
  970 ha(r1,j)=r6
  980 continue
  990 continue
1060  continue
      aflag(5)=aflag(7)/aflag(6)
      iflag(1)=-2
 1110 z=zz
      return
      end subroutine














      subroutine y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)
      implicit double precision(a-b,g,p,t-y),integer (c,f,h-n,r-s,z)
      double precision a(nn), pivot(n), b(n)
      integer snr(nn), ha(iha,11), iflag(10)
      ifail=0
      if(iflag(1).eq.-2) goto 1000
      ifail=1
      goto 1110
1000  mode=iflag(4)
      ipiv=iflag(3)
      n8=n+1
      n7=n-1
      state=iflag(5)
!
!  solve the system with lower triangular matrix  l  (if the
!  lu-factorization is available).
!
      if(state.ne.3) goto 1051
      if(ipiv.eq.0) goto 1020
      do 1010 i=1,n7
      l1=ha(i,7)
      t=b(l1)
      b(l1)=b(i)
      b(i)=t
1010  continue
1020  continue
      do 1050 i=1,n
      rr1=ha(i,1)
      rr2=ha(i,2)-1
      if(rr1.gt.rr2) goto 1040
      do 1030 j=rr1,rr2
      l1=snr(j)
 1030 b(i)=b(i)-a(j)*b(l1)
 1040 continue
 1050 continue
!
!  solve the system with upper triagular matrix.
!
 1051 continue
      do 1090 i=1,n
      r1=n8-i
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      if(rr2.lt.rr1)   goto 1080
      do 1070 j=rr1,rr2
      r2=snr(j)
 1070 b(r1)=b(r1)-a(j)*b(r2)
 1080 continue
 1090 b(r1)=b(r1)/pivot(r1)
!
! if interchanges were used during the  elimination then a reordering in
! lution vector is made.
!
      if(ipiv.eq.0) goto 1110
      do 1100 i=1,n7
      r1=n-i
      r2=ha(r1,8)
      t=b(r2)
      b(r2)=b(r1)
 1100 b(r1)=t
 1110 return
      end subroutine

end module sparse_mat
