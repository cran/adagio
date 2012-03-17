*** ********************************************************************
      subroutine risort (a, m, n)
c-----------------------------------------------------------------------
c     the shell sorting procedure is used to reorder the elements of a
c     so that a(i).le.a(i+1) for i=1,...,n-1. the same permutations are
c     performed on m that are performed on a. it is assumed that n.ge.1.
c-----------------------------------------------------------------------
      real a(n)
      integer m(n), t
      integer k(10)
c------------------------
      data k(1)/1/, k(2)/4/, k(3)/13/, k(4)/40/, k(5)/121/, k(6)/364/,
     1     k(7)/1093/, k(8)/3280/, k(9)/9841/, k(10)/29524/
c------------------------
c
c             selection of the increments k(i) = (3**i-1)/2
c
      if (n .lt. 2) return
      imax = 1
      do 10 i = 3,10
         if (n .le. k(i)) go to 20
         imax = imax + 1
   10 continue
c
c            stepping through the increments k(imax),...,k(1)
c
   20 i = imax
      do 40 ii = 1,imax
         ki = k(i)
c
c             sorting elements that are ki positions apart
c               so that a(j).le.a(j+ki) for j=1,...,n-ki
c
         jmax = n - ki
         do 32 j = 1,jmax
            l = j
            ll = j + ki
            s = a(ll)
            t = m(ll)
   30          if (s .ge. a(l)) go to 31
               a(ll) = a(l)
               m(ll) = m(l)
               ll = l
               l = l - ki
               if (l .gt. 0) go to 30
   31       a(ll) = s
            m(ll) = t
   32    continue
c
   40 i = i - 1
      return
      end

*** ********************************************************************
      subroutine sigma1 (n,m,p,w,q,i,b,kub,ub,np1,n5,lx,lr,
     *                   bs,ps,ws,xs,iwk)
c
c subroutine to compute an upper bound  ub  on the best
c final solution which can be obtained from the current
c solution.
c
      integer p(n),w(n),q(m),b(np1),ub,iwk(n5)
      integer lx(n),bs(n),ps(np1),ws(np1),xs(n)
      integer qs,sb
c
      ns = 0
      qs = 0
      do 10 j=i,m
        qs = qs + q(j)
   10 continue
      sb = 0
      do 20 j=1,n
        lx(j) = 0
        if ( b(j) .eq. 0 ) go to 20
        ns = ns + 1
        bs(ns) = j
        ps(ns) = p(j)
        ws(ns) = w(j)
        sb = sb + w(j)
   20 continue
      if ( sb .gt. qs ) go to 40
      lr = qs - sb
      ub = 0
      if ( ns .eq. 0 ) return
      do 30 j=1,ns
        ub = ub + ps(j)
        xs(j) = 1
   30 continue
      go to 50
   40 call sknp (ns,qs,kub,ub,n,np1,n5,ps,ws,xs,iwk)
      lr = qs
   50 do 60 j=1,ns
        jj = bs(j)
        lx(jj) = xs(j)
   60 continue
      return
      end

*** ********************************************************************
      subroutine pi1 (n,m,p,w,q,i,b,bb,kub,bl,lb,pbl,v,xl,
     *                np1,n5,bs,ps,ws,xs,iwk)
c
c subroutine to compute a feasible solution to the current
c problem. the solution is stored in array  xl , the
c corresponding value in  lb .
c
      integer bb(m,n),bl(m,np1),xl(m,n),iwk(n5)
      integer p(n),w(n),q(m),b(np1),pbl(m),v(m)
      integer bs(n),ps(np1),ws(np1),xs(n)
      integer pb,qs,sb,u
c
c step 1
c
      u = 0
      do 10 j=1,n
        if ( b(j) .eq. 0 ) go to 10
        u = u + 1
        bs(u) = j
   10 continue
      do 20 j=i,m
        pbl(j) = 0
        v(j) = 0
   20 continue
      lb = 0
      ikub = kub
      if ( u .eq. 0 ) return
      ns = 0
      sb = 0
      do 30 j=1,u
        jj = bs(j)
        if ( bb(i,jj) .ne. 0 ) go to 30
        if ( w(jj) .gt. q(i) ) go to 30
        ns = ns + 1
        sb = sb + w(jj)
        bl(i,ns) = jj
        ps(ns) = p(jj)
        ws(ns) = w(jj)
   30 continue
      ii = i
c
c step 2
c
   40 pbl(ii) = ns
      if ( sb .gt. q(ii) ) go to 60
      pb = 0
      if ( ns .eq. 0 ) go to 80
      do 50 j=1,ns
        pb = pb + ps(j)
        xl(ii,j) = 1
   50 continue
      go to 80
   60 qs = q(ii)
      kub = 0
      if ( ii .eq. m ) kub = ikub
      call sknp (ns,qs,kub,pb,n,np1,n5,ps,ws,xs,iwk)
      do 70 j=1,ns
        xl(ii,j) = xs(j)
   70 continue
   80 lb = lb + pb
      ikub = ikub - pb
      v(ii) = pb
      bl(ii,ns+1) = n + 1
c
c step 3
c
      if ( ii .eq. m ) return
      jb = 1
      jbs = 0
      do 100 j=1,u
        if ( bs(j) .lt. bl(ii,jb) ) go to 90
        jb = jb + 1
        if ( xl(ii,jb-1) .eq. 1 ) go to 100
   90   jbs = jbs + 1
        bs(jbs) = bs(j)
  100 continue
      u = jbs
      if ( u .eq. 0 ) return
      ns = 0
      sb = 0
      ii = ii + 1
      do 110 j=1,u
        jj = bs(j)
        if( w(jj) .gt. q(ii) ) go to 110
        ns = ns + 1
        sb = sb + w(jj)
        bl(ii,ns) = jj
        ps(ns) = p(jj)
        ws(ns) =  w(jj)
  110 continue
      go to 40
      end

*** ********************************************************************
      subroutine parc (i,ii,ub,iflag,vb,lub,lj,li,f,bb,q,b,n,m,np1,
     *                 lx,lxi,lr,lri,lubi)
c
c subroutine for parametric computation of the upper bounds.
c
      integer f(m),bb(m,n),q(m),b(np1),ub,vb,r,s
      integer lx(n),lxi(n)
c
      iflag = 0
      if ( b(lj) .ne. 0 ) go to 60
      i1 = i - 1
      if ( i1 .lt. li ) go to 20
      iq = 0
      do 10 r=li,i1
        iq = iq + q(r)
   10 continue
      if ( iq .gt. lr ) return
   20 r = ii
      s = f(r)
   30 if ( s .ne. (-1) ) go to 40
      r = r - 1
      s = f(r)
      go to 30
   40 if ( lx(s) .eq. 0 ) return
      if ( s .eq. lj ) go to 50
      s = bb(r,s)
      go to 30
   50 ub = lub - vb
      iflag = 1
      return
   60 i1 = i - 1
      if ( i1 .lt. 1 ) go to 80
      iq = 0
      do 70 r=1,i1
        iq = iq + q(r)
   70 continue
      if ( iq .gt. lri ) return
   80 do 90 j=1,n
        if ( b(j) .eq. 1 ) go to 90
        if ( lxi(j) .eq. 0 ) return
   90 continue
      ub = lubi - vb
      iflag = 1
      return
      end

*** ********************************************************************
      subroutine sknp (ns,qs,kub,vs,n,np1,n5,ps,ws,xs,iwk)
c
c subroutine to solve the 0-1 single knapsack problem
c
c maximize    vs = ps(1)*xs(1) + ... + ps(ns)*xs(ns)
c subject to       ws(1)*xs(1) + ... + ws(ns)*xs(ns) .le. qs
c                  xs(j) = 0 or 1   for  j=1,...,ns
c                  vs .gt. kub
c
c this subroutine is a modified version of subroutine kp01
c which appeared in  computing 21, 81-86(1978).
c
      integer qs, vs
      integer ps(np1), ws(np1), xs(n), iwk(n5)
c
      i1 = 1
      i2 = i1 + n
      i3 = i2 + n
      i4 = i3 + n
      i5 = i4 + n
      call sknp1 (ns,qs,kub,vs,n,np1,ps,ws,xs,iwk(i1),iwk(i2),
     *            iwk(i3),iwk(i4),iwk(i5))
      return
      end

*** ********************************************************************
      subroutine sknp1 (ns,qs,kub,vs,n,np1,ps,ws,xs,d,min,
     *                  pbar,wbar,zbar)
c
c subroutine to solve the 0-1 single knapsack problem
c
c maximize    vs = ps(1)*xs(1) + ... + ps(ns)*xs(ns)
c subject to       ws(1)*xs(1) + ... + ws(ns)*xs(ns) .le. qs
c                  xs(j) = 0 or 1   for  j=1,...,ns
c                  vs .gt. kub
c
c this subroutine is a modified version of subroutine kp01
c which appeared in  computing 21, 81-86(1978).
c
      integer qs,vs,diff,pr,r,t
      integer ps(np1),ws(np1),xs(n)
      integer d(n),min(n),pbar(n),wbar(n),zbar(n)
c
      vs = kub
      ip = 0
      ms = qs
      do 10 l=1,ns
        ll = l
        if ( ws(l) .gt. ms ) go to 20
        ip = ip + ps(l)
        ms = ms - ws(l)
   10 continue
   20 ll = ll - 1
      if ( ms .eq. 0 ) go to 50
      ps(ns+1) = 0
      ws(ns+1) = qs + 1
      lim = ip + (ms*ps(ll+2))/ws(ll+2)
      a = ip + ps(ll+1)
      b = (ws(ll+1) - ms)*ps(ll)
      c = ws(ll)
      lim1 = a - b/c
      if ( lim1 .gt. lim ) lim = lim1
      if ( lim .le. vs ) return
      mink = qs + 1
      min(ns) = mink
      do 30 j=2,ns
        kk = ns + 2 - j
        if ( ws(kk) .lt. mink ) mink = ws(kk)
        min(kk-1) = mink
   30 continue
      do 40 j=1,ns
        d(j) = 0
   40 continue
      pr = 0
      lold = ns
      ii = 1
      go to 170
   50 if ( vs .ge. ip ) return
      vs = ip
      do 60 j=1,ll
        xs(j) = 1
   60 continue
      nn = ll + 1
      do 70 j=nn,ns
        xs(j) = 0
   70 continue
      qs = 0
      return
   80 if ( ws(ii) .le. qs ) go to 90
      ii1 = ii + 1
      if ( vs .ge. (qs*ps(ii1))/ws(ii1) + pr ) go to 280
      ii = ii1
      go to 80
   90 ip = pbar(ii)
      ms = qs - wbar(ii)
      in = zbar(ii)
      ll = ns
      if ( in .gt. ns) go to 110
      do 100 l=in,ns
        ll = l
        if ( ws(l) .gt. ms ) go to 160
        ip = ip + ps(l)
        ms = ms - ws(l)
  100 continue
  110 if ( vs .ge. ip + pr ) go to 280
      vs = ip + pr
      mfirst = ms
      nn = ii - 1
      do 120 j=1,nn
        xs(j) = d(j)
  120 continue
      do 130 j=ii,ll
        xs(j) = 1
  130 continue
      if ( ll .eq. ns ) go to 150
      nn = ll + 1
      do 140 j=nn,ns
        xs(j) = 0
  140 continue
  150 if ( vs .ne. lim ) go to 280
      qs = mfirst
      return
  160 l = ll
      ll = ll - 1
      if ( ms .eq. 0 ) go to 110
      if ( vs .ge. pr + ip + (ms*ps(l))/ws(l) ) go to 280
  170 wbar(ii) = qs - ms
      pbar(ii) = ip
      zbar(ii) = ll + 1
      d(ii) = 1
      nn = ll - 1
      if ( nn .lt. ii ) go to 190
      do 180 j=ii,nn
        wbar(j+1) = wbar(j) - ws(j)
        pbar(j+1) = pbar(j) - ps(j)
        zbar(j+1) = ll + 1
        d(j+1) = 1
  180 continue
  190 j1 = ll + 1
      do 200 j=j1,lold
        wbar(j) = 0
        pbar(j) = 0
        zbar(j) = j
  200 continue
      lold = ll
      qs = ms
      pr = pr + ip
      if ( ll - (ns - 2) ) 240, 220, 210
  210 ii = ns
      go to 250
  220 if ( qs .lt. ws(ns) ) go to 230
      qs = qs - ws(ns)
      pr = pr + ps(ns)
      d(ns) = 1
  230 ii = ns - 1
      go to 250
  240 ii = ll + 2
      if ( qs .ge. min(ii-1) ) go to 80
  250 if ( vs .ge. pr ) go to 270
      vs = pr
      do 260 j=1,ns
        xs(j) = d(j)
  260 continue
      mfirst = qs
      if ( vs .eq. lim ) return
  270 if ( d(ns) .eq. 0 ) go to 280
      d(ns) = 0
      qs = qs + ws(ns)
      pr = pr - ps(ns)
  280 nn = ii - 1
      if ( nn .eq. 0 ) go to 300
      do 290 j=1,nn
        kk = ii - j
        if ( d(kk) .eq. 1 ) go to 310
  290 continue
  300 qs = mfirst
      return
  310 r = qs
      qs = qs + ws(kk)
      pr = pr - ps(kk)
      d(kk) = 0
      if ( r .lt. min(kk) ) go to 320
      ii = kk + 1
      go to 80
  320 nn = kk + 1
      ii = kk
  330 if ( vs .ge. pr + (qs*ps(nn))/ws(nn) ) go to 280
      diff = ws(nn) - ws(kk)
      if ( diff ) 390, 340, 350
  340 nn = nn + 1
      go to 330
  350 if ( diff .gt. r ) go to 340
      if ( vs .ge. pr + ps(nn) ) go to 340
      vs = pr + ps(nn)
      do 360 j=1,kk
        xs(j) = d(j)
  360 continue
      jj = kk + 1
      do 370 j=jj,ns
        xs(j) = 0
  370 continue
      xs(nn) = 1
      mfirst = qs - ws(nn)
      if ( vs .ne. lim ) go to 380
      qs = mfirst
      return
  380 r = r - diff
      kk = nn
      nn = nn + 1
      go to 330
  390 t = r - diff
      if ( t .lt. min(nn) ) go to 340
      n1 = nn + 1
      if ( vs .ge. pr + ps(nn) + (t*ps(n1))/ws(n1) ) go to 280
      qs = qs - ws(nn)
      pr = pr + ps(nn)
      d(nn) = 1
      ii = nn + 1
      wbar(nn) = ws(nn)
      pbar(nn) = ps(nn)
      zbar(nn) = ii
      do 400 j=n1,lold
        wbar(j) = 0
        pbar(j) = 0
        zbar(j) = j
  400 continue
      lold = nn
      go to 80
      end
