*** ********************************************************************
***     S u b r o u t i n e  M K P
*** ********************************************************************

      subroutine mkp (n,m,p,w,k,bck,xstar,vstar,wk,iwk,num)
*   ____________________________________________________________________
*
*   Subroutine to solve a 0-1 multiple knapsack problem of n
*   items (with  n .ge. 2) and  m  knapsacks (with  m .ge. 1 ,
*   i.e. also suitable for a 0-1 single knapsack problem).
*   The problem to be solved is
*
*   maximize  vstar = p(1)*(x(1,1) + ... + x(m,1)) + ...
*                     ... + p(n)*(x(1,n) + ... + x(m,n))
*   subject to
*   w(1)*x(i,1) + ... + w(n)*x(i,n) .le. k(i)   for  i=1,...,m
*   x(1,j) + ... + x(m,j) .le. 1   for  j=1,...,n
*   x(i,j) = 0 or 1   for  i=1,...,m ,  j=1,...,n ,
*
*   where all p(j), w(j), and k(i) are positive integers.
*   Before mkp is called, array k must be sorted so that
*   k(1) .le. k(2) .le. ... .le. k(m) .
*
*   Meaning of the input parameters ...
*
*   n    = number of items.
*   m    = number of knapsacks.
*   p(j) = profit of item  j  (j=1,...,n) .
*   w(j) = weight of item  j  (j=1,...,n) .
*   k(i) = capacity of knapsack  i  (i=1,...,m) .
*   bck  = -1  if exact solution is required.
*        = maximum number of backtrackings to be performed, if
*          heuristic solution is required.
*   wk   = real work space of dimension n.
*   iwk  = work space of dimension .ge. 5*m + 14*n + 4*m*n + 3
*   num  = dimension of iwk
*
*   Meaning of the output parameters ...
*
*   xstar(j) = 0  if item  j  is not in the optimal solution
*                 (i.e. if  x(i,j) = 0  for all  i ).
*            = knapsack where item  j  is inserted, otherwise
*              (i.e. if  x(xstar(j),j) = 1 ).
*   vstar    = value of the optimal solution if  vstar .gt. 0.
*            = error condition (infeasibility or triviality)
*              in the input data if  vstar .lt. 0 .
*              = -1  if  n .lt. 2  or  m .lt. 1 .
*              = -2  if some  p(j) ,  w(j)  or  k(i) are not
*                    positive.
*              = -3  if a knapsack cannot contain any item.
*              = -4  if an item cannot fit into any knapsack.
*              = -5  if knapsack  m  contains all the items.
*              = -7  if array  k is not correctly sorted.
*              = -8  if num .lt. 5*m + 14*n + 4*m*n + 3.
*  
*              (in all the above cases array  xstar is not
*              defined).
*
*   All the parameters except wk are of integer type. when mkp
*   terminates, all the input parameters are unchanged except
*   bck, which gives the number of backtrackings performed.
*   ____________________________________________________________________

      integer p(n),w(n),k(m),xstar(n),bck,vstar,iwk(num)
      real wk(n)
      integer bb, bl, x, xl
      integer b, ubb
      integer f, pbl, q, v
      integer bs, ps, ws, xs
C
C                    check the input data
C
      if (m .lt. 1 .or. n .lt. 2) go to 100
      mn = m*n
      if (num .lt. 5*m + 14*n + 4*mn + 3) go to 160

      if (p(1) .le. 0 .or. w(1) .le. 0) go to 110
      ap = p(1)
      aw = w(1)
      wk(1) = -ap/aw
      maxw = w(1)
      minw = w(1)
      isumw = w(1)
      do 10 j = 2,n
         if (p(j) .le. 0 .or. w(j) .le. 0) go to 110
         ap = p(j)
         aw = w(j)
         wk(j) = -ap/aw
         if (w(j) .gt. maxw) maxw = w(j)
         if (w(j) .lt. minw) minw = w(j)
         isumw = isumw + w(j)
   10 continue

      if (k(1) .le. 0) go to 110
      if (m .eq. 1) go to 30
      do 20 i = 2,m
         if (k(i) .le. 0) go to 110
         if (k(i) .lt. k(i-1)) go to 150
   20 continue

   30 if (minw .gt. k(1)) go to 120
      if (maxw .gt. k(m)) go to 130
      if (isumw .le. k(m)) go to 140
      vstar = 0
C
C             reorder the arrays p and w so that
C                p(j)/w(j) .ge. p(j+1)/w(j+1)
C
      n5 = 5*n
      do 40 j = 1,n
         jj = n5 + j
         iwk(jj) = j
   40 continue
      call risort (wk, iwk(n5 + 1), n)

      do 50 j = 1,n
         iwk(j) = p(j)
         jn = j + n
         iwk(jn) = w(j)
   50 continue

      do 60 j = 1,n
         jj = n5 + j
         l = iwk(jj)
         p(j) = iwk(l)
         npl = n + l
         w(j) = iwk(npl)
   60 continue
C
C                partition the work space iwk
C
      lx =  jj + 1
      lxi = lx + n
      bs = lxi + n
      xs =  bs + n
      ubb = xs + n

      np1 = n + 1
      b = ubb + n
      ps = b  + np1
      ws = ps + np1

      f = ws + np1
      pbl = f + m
      q = pbl + m
      v = q + m

      bb = v + m
      x = bb + mn
      xl = x + mn

      bl = xl + mn
C
C                     solve the problem
C
      call mkp1 (n, m, p, w, k, bck, xstar, vstar, np1, n5,
     +           iwk(bb), iwk(bl), iwk(x), iwk(xl),
     +           iwk(b), iwk(ubb), iwk(lx), iwk(lxi),
     +           iwk(f), iwk(pbl), iwk(q), iwk(v),
     +           iwk(bs), iwk(ps), iwk(ws), iwk(xs), iwk(1))
C
C           restore the initial ordering to p and w,
C                and reorder xstar accordingly
C
      do 70 j = 1,n
         iwk(j) = p(j)
         jn = j + n
         iwk(jn) = w(j)
         jnn = jn + n
         iwk(jnn) = xstar(j)
   70 continue

      do 80 j = 1,n
         jj = n5 + j
         l = iwk(jj)
         p(l) = iwk(j)
         jn = j + n
         w(l) = iwk(jn)
         jnn = jn + n
         xstar(l) = iwk(jnn)
   80 continue
      return
C
C                        error return
C
  100 vstar = -1
      return
  110 vstar = -2
      return
  120 vstar = -3
      return
  130 vstar = -4
      return
  140 vstar = -5
      return
  150 vstar = -7
      return
  160 vstar = -8
      return
      end


*** ********************************************************************
***     S u b r o u t i n e  M K P 1
*** ********************************************************************

      subroutine mkp1 (n, m, p, w, k, bck, xstar, vstar, np1, n5,
     +                 bb, bl, x, xl, b, ubb, lx, lxi,
     +                 f, pbl, q, v, bs, ps, ws, xs, iwk)

*   --------------------------------------------------------------------
*   meaning of the main internal variables and arrays ...
*  
*   i       = knapsack currently considered.
*   lb      = lower bound on the optimal solution.
*   ub      = upper bound on the optimal solution.
*   vb      = value of the current solution.
*   x(i,j)  = 1  if item  j  is inserted in knapsack  i  in
*                the current solution.
*           = 0  otherwise.
*   f(i)    = pointer to the last item inserted in knapsack  i
*             ( = -1  if knapsack  i  is empty).
*   bb(i,j) = pointer to the item inserted in knapsack  i
*             just before item  j ( = -1  if  j  is the first
*             item inserted in knapsack  i ).
*   q(i)    = current available capacity of knapsack  i .
*   b(j)    = 1  if item  j  is not inserted in any knapsack.
*           = 0  if item  j  is inserted in a knapsack.
*   pbl(i)  = number of the items which can be inserted in
*             knapsack  i .
*   bl(i,s) = pointer to the  s-th  item which can be inserted
*             in knapsack  i .
*   xl(i,j) = 1  if item  j  was inserted in knapsack  i  in
*                the last execution of subroutine pi1.
*           = 0  otherwise.
*   iwk       work space for the subroutine sknp.
*   --------------------------------------------------------------------

      integer p(n), w(n), k(m), xstar(n), bck, vstar
      integer bb(m,n), bl(m,np1), x(m,n), xl(m,n)
      integer b(np1), ubb(n), lx(n), lxi(n)
      integer f(m), pbl(m), q(m), v(m)
      integer bs(n), ps(np1), ws(np1), xs(n), iwk(n5)
      integer s, u, ub, vb

      if (m .eq. 1) go to 250
C
C                   step 1 (initialization)
C
      jbck = bck
      bck = 0
      kub = 0
      n1 = n + 1
      b(n1) = 1
      m1 = m - 1
      do 40 j=1,n
        b(j) = 1
        do 30 i=1,m
          x(i,j) = 0
          bb(i,j) = 0
   30   continue
   40 continue
      do 50 i=1,m1
        q(i) = k(i)
        f(i) = -1
   50 continue
      q(m) = k(m)
      vstar = 0
      vb = 0
      i = 1
      call sigma1 (n,m,p,w,k,1,b,kub,ub,np1,n5,lx,lr,
     +             bs,ps,ws,xs,iwk)
      do 60 j=1,n
        lxi(j) = lx(j)
   60 continue
      lri = lr
      lubi = ub
      iflag = 0
C
C                   step 2 (heuristic)
C
   70 kub = vstar - vb
      call pi1 (n,m,p,w,q,i,b,bb,kub,bl,lb,pbl,v,xl,
     +          np1,n5,bs,ps,ws,xs,iwk)
      if ( lb + vb .le. vstar ) go to 140
      vstar = lb + vb
      do 90 j=1,n
        xstar(j) = 0
        do 80 s=1,i
          if ( x(s,j) .eq. 0 ) go to 80
          xstar(j) = s
          go to 90
   80   continue
   90 continue
      ip = pbl(i)
      if ( ip .eq. 0 ) go to 110
      do 100 j=1,ip
        jj = bl(i,j)
        if ( xl(i,j) .eq. 1 ) xstar(jj) = i
  100 continue
  110 i1 = i + 1
      do 130 ii=i1,m
        ip = pbl(ii)
        if ( ip .eq. 0 ) go to 130
        do 120 j=1,ip
          jj = bl(ii,j)
          if ( xl(ii,j) .eq. 1 ) xstar(jj) = ii
  120   continue
  130 continue
      if ( ub .eq. lb ) go to 200
C
C                   step 3 (updating)
C
  140 if ( v(i) .eq. 0 ) go to 180
      iuv = ub + vb
      u = pbl(i)
      ibv = 0
      do 170 s=1,u
        if ( xl(i,s) .eq. 0 ) go to 170
        j = bl(i,s)
        x(i,j) = 1
        q(i) = q(i) - w(j)
        vb = vb + p(j)
        b(j) = 0
        bb(i,j) = f(i)
        ubb(j) = iuv
        if ( iflag .eq. 1 ) go to 150
        lub = iuv
        lj = j
        li = i
  150   f(i) = j
        ibv = ibv + p(j)
        if ( ibv .eq. v(i) ) go to 180
        call parc (i,i,ub,iflag,vb,lub,lj,li,f,bb,q,b,n,m,np1,
     +             lx,lxi,lr,lri,lubi)
        if ( iflag .eq. 1 ) go to 160
        kub = vstar - vb
        call sigma1 (n,m,p,w,q,i,b,kub,ub,np1,n5,lx,lr,
     +               bs,ps,ws,xs,iwk)
        lj = n1
  160   iuv = ub + vb
        if ( iuv .le. vstar ) go to 200
  170 continue
  180 if ( i .eq. m - 1 ) go to 200
      ip1 = i + 1
      call parc (ip1,i,ub,iflag,vb,lub,lj,li,f,bb,q,b,n,m,np1,
     +           lx,lxi,lr,lri,lubi)
      if ( iflag .eq. 1 ) go to 190
      kub = vstar - vb
      call sigma1 (n,m,p,w,q,ip1,b,kub,ub,np1,n5,lx,lr,
     +             bs,ps,ws,xs,iwk)
      lj = n1
  190 if ( ub + vb .le. vstar ) go to 200
      i = i + 1
      go to 140
C
C                   step 4 (backtracking)
C
  200 if ( i .gt. 0 ) go to 210
      bck = bck - 1
      return
  210 if ( bck .eq. jbck ) return
      bck = bck + 1
      if ( f(i) .ne. (-1) ) go to 230
      do 220 j=1,n
        bb(i,j) = 0
  220 continue
      i = i - 1
      go to 200
  230 j = f(i)
      x(i,j) = 0
      b(j) = 1
      vb = vb - p(j)
      q(i) = q(i) + w(j)
      do 240 s=1,n
        if ( bb(i,s) .eq. j ) bb(i,s) = 0
  240 continue
      f(i) = bb(i,j)
      if ( ubb(j) .le. vstar ) go to 200
      ub = ubb(j) - vb
      iflag = 1
      go to 70
C
C           particular case (0-1 single knapsack problem)
C
  250 k1 = k(1)
      do 260 j=1,n
        ps(j) = p(j)
        ws(j) = w(j)
  260 continue
      call sknp (n,k1,0,vstar,n,np1,n5,ps,ws,xs,iwk)
      do 270 j=1,n
        xstar(j) = xs(j)
  270 continue
      bck = 0
      return
      end
