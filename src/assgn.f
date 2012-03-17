*** ********************************************************************
      subroutine assgn (n,a,c,t,iwk,ierr)
C     -------------------
C     solution of the assignment problem
C     -------------------
      integer a(n,*), c(n), t, iwk(*)

      i1 = n + 1
      i2 = i1 + n
      i3 = i2 + n
      i4 = i3 + n + 1
      i5 = i4 + n
      i6 = i5 + n
      call assgn1(n,a,c,t,iwk(1),iwk(i1),iwk(i2),iwk(i3),iwk(1),
     +            iwk(i3),iwk(i4),iwk(i5),iwk(i6),ierr)
      return
      end

*** ********************************************************************
      subroutine assgn1(n,a,c,t,ch,lc,lr,lz,nz,rh,slc,slr,
     +                  u,ierr)
      integer a(n,*), c(n), ch(n), lc(n), lr(n), lz(n),
     +        nz(n), rh(*), slc(n), slr(n), u(*)
      integer h, q, r, s, t
C    ___________________________________________________________________
C
C    This subroutine solves the square assignment problem
C    The meaning of the input parameters is
C    n = number of rows and columns of the cost matrix
C    a(i,j) = element in row i and column j of the cost matrix
C    ( at the end of computation the elements of a are changed)
C    the meaning of the output parameters is
C    c(j) = row assigned to column j (j=1,n)
C    t = cost of the optimal assignment
C    All parameters are integer
C    The meaning of the local variables is
C    a(i,j) = element of the cost matrix if a(i,j) is positive,
C             column of the unassigned zero following in row i
C             (i=1,n) the unassigned zero of column j (j=1,n)
C             if a(i,j) is not positive
C    a(i,n+1) = column of the first unassigned zero of row i
C               (i=1,n)
C    ch(i) = column of the next unexplored and unassigned zero
C            of row i (i=1,n)
C    lc(j) = label of column j (j=1,n)
C    lr(i) = label of row i (i=1,n)
C    lz(i) = column of the last unassigned zero of row i(i=1,n)
C    nz(i) = column of the next unassigned zero of row i(i=1,n)
C    rh(i) = unexplored row following the unexplored row i
C            (i=1,n)
C    rh(n+1) = first unexplored row
C    slc(k) = k-th element contained in the set of the labelled
C             columns
C    slr(k) = k-th element contained in the set of the labelled
C             rows
C    u(i) = unassigned row following the unassigned row i
C           (i=1,n)
C    u(n+1) = first unassigned row
C    ierr = 0 if the routine terminates successfully. otherwise
C             ierr = 1
C   
C    The vectors c,ch,lc,lr,lz,nz,slc,slr must be dimensioned
C    at least at (n), the vectors rh,u at  least at (n+1),
C    and the matrix a at least at (n,n+1). To save storage
C    lz and rh may use the same storage area, and nz and ch
C    may use the same storage area.
C    ___________________________________________________________________

C initialization
      maxnum = 2147483647
      ierr = 0
      np1 = n+1
      do 10 j=1,n
        c(j) = 0
        lz(j) = 0
        nz(j) = 0
        u(j) = 0
   10 continue
      u(np1) = 0
      t = 0
C reduction of the initial cost matrix
      do 40 j=1,n
        s = a(1,j)
        do 15 l=2,n
          if ( a(l,j) .lt. s ) s = a(l,j)
   15   continue
        if (s) 20,40,30
   20   mm = maxnum + s
        if (t .lt. -mm) go to 400
        t = t + s
        do 25 i = 1,n
          if (a(i,j) .gt. mm) go to 400
          a(i,j) = a(i,j) - s
   25   continue
        go to 40
   30   mm = maxnum - s
        if (t .gt. mm) go to 400
        t = t + s
        do 35 i = 1,n
          a(i,j) = a(i,j) - s
   35   continue
   40 continue
      do 70 i=1,n
        q = a(i,1)
        do 50 l=2,n
          if ( a(i,l) .lt. q ) q = a(i,l)
   50   continue
        mm = maxnum - q
        if (t .gt. mm) go to 400
        t = t + q
        l = np1
        do 60 j=1,n
          a(i,j) = a(i,j)-q
          if ( a(i,j) .ne. 0 ) go to 60
          a(i,l) = -j
          l = j
   60   continue
   70 continue
C choice of the initial solution
      k = np1
      do 140 i=1,n
        lj = np1
        j = -a(i,np1)
   80   if ( c(j) .eq. 0 ) go to 130
        lj = j
        j = -a(i,j)
        if ( j .ne. 0 ) go to 80
        lj = np1
        j = -a(i,np1)
   90   r = c(j)
        lm = lz(r)
        m = nz(r)
  100   if ( m .eq. 0 ) go to 110
        if ( c(m) .eq. 0 ) go to 120
        lm = m
        m = -a(r,m)
        go to 100
  110   lj = j
        j = -a(i,j)
        if ( j .ne. 0 ) go to 90
        u(k) = i
        k = i
        go to 140
  120   nz(r) = -a(r,m)
        lz(r) = j
        a(r,lm) = -j
        a(r,j) = a(r,m)
        a(r,m) = 0
        c(m) = r
  130   c(j) = i
        a(i,lj) = a(i,j)
        nz(i) = -a(i,j)
        lz(i) = lj
        a(i,j) = 0
  140 continue
C research of a new assignment
  150 if ( u(np1) .eq. 0 ) return
      do 160 i=1,n
        ch(i) = 0
        lc(i) = 0
        lr(i) = 0
        rh(i) = 0
  160 continue
      rh(np1) = -1
      kslc = 0
      kslr = 1
      r = u(np1)
      lr(r) = -1
      slr(1) = r
      if ( a(r,np1) .eq. 0 ) go to 220
  170 l = -a(r,np1)
      if ( a(r,l) .eq. 0 ) go to 180
      if ( rh(r) .ne. 0 ) go to 180
      rh(r) = rh(np1)
      ch(r) = -a(r,l)
      rh(np1) = r
  180 if ( lc(l) .eq. 0 ) go to 200
      if ( rh(r) .eq. 0 ) go to 210
  190 l = ch(r)
      ch(r) = -a(r,l)
      if ( a(r,l) .ne. 0 ) go to 180
      rh(np1) = rh(r)
      rh(r) = 0
      go to 180
  200 lc(l) = r
      if ( c(l) .eq. 0 ) go to 360
      kslc = kslc+1
      slc(kslc) = l
      r = c(l)
      lr(r) = l
      kslr = kslr+1
      slr(kslr) = r
      if ( a(r,np1) .ne. 0 ) go to 170
  210 continue
      if ( rh(np1) .gt. 0 ) go to 350
C reduction of the current cost matrix
  220 h = maxnum
      do 240 j=1,n
        if ( lc(j) .ne. 0 ) go to 240
        do 230 k=1,kslr
          i = slr(k)
          if ( a(i,j) .lt. h ) h = a(i,j)
  230   continue
  240 continue
      mm = maxnum - h
      if (mm .eq. 0 .or. t .gt. mm) go to 400
      t = t + h
      do 290 j=1,n
        if ( lc(j) .ne. 0 ) go to 290
        do 280 k=1,kslr
          i = slr(k)
          a(i,j) = a(i,j)-h
          if ( a(i,j) .ne. 0 ) go to 280
          if ( rh(i) .ne. 0 ) go to 250
          rh(i) = rh(np1)
          ch(i) = j
          rh(np1) = i
  250     l = np1
  260     nl = -a(i,l)
          if ( nl .eq. 0 ) go to 270
          l = nl
          go to 260
  270     a(i,l) = -j
  280   continue
  290 continue
      if ( kslc .eq. 0 ) go to 350
      do 340 i=1,n
        if ( lr(i) .ne. 0 ) go to 340
        do 330 k=1,kslc
          j = slc(k)
          if ( a(i,j) .gt. 0 ) go to 320
          l = np1
  300     nl = - a(i,l)
          if ( nl .eq. j ) go to 310
          l = nl
          go to 300
  310     a(i,l) = a(i,j)
          a(i,j) = h
          go to 330
  320     mm = maxnum - h
          if (a(i,j) .gt. mm) go to 400
          a(i,j) = a(i,j) + h
  330   continue
  340 continue
  350 r = rh(np1)
      go to 190
C assignment of a new row
  360 c(l) = r
      m = np1
  370 nm = -a(r,m)
      if ( nm .eq. l ) go to 380
      m = nm
      go to 370
  380 a(r,m) = a(r,l)
      a(r,l) = 0
      if ( lr(r) .lt. 0 ) go to 390
      l = lr(r)
      a(r,l) = a(r,np1)
      a(r,np1) = -l
      r = lc(l)
      go to 360
  390 u(np1) = u(r)
      u(r) = 0
      go to 150
C error return - integer overflow occurs
  400 ierr = 1
      return
      end
