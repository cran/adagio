C Maximal Sum Subarray Problem
      subroutine maxsubf(x, n, s, i1, i2)
      integer n, i1, i2, j1, j2
      double precision x(n)
      double precision s, ss

      ss = 0
      j1 = 1
      j2 = 1

      do 3000 i = 1,n 
      if (ss .gt. -x(i)) then
        ss = ss + x(i)
        j2 = i
        if (ss .gt. s) then
          s  = ss
          i1 = j1
          i2 = j2
        endif
      else
        ss = 0
        j1 = i+1
        j2 = i+1
      endif
 3000 continue

      return
      end


C Maximal Sum Subrectangle Problem
      subroutine maxsub2f(a, s, n, m, fmax, mind, aa, b)
      integer n, m, mind(4)
      integer mi1, mi2
      double precision a(n, m), aa(m), b, s(n+1, m)
      double precision fsum, fmax

      do 3000 i = 1, m 
        s(1, i) = 0
 3000 continue

      do 3002 i = 2, n+1 
        b = a(i-1, 1)
        s(i, 1) = s(i-1, 1) + b
        do 3004 j = 2, m 
          b = b + a(i-1, j)
          s(i, j) = s(i-1, j) + b
 3004   continue
 3002 continue

      fsum = 0.0
      mi1 = 0
      mi2 = 0

      do 3006 i = 2, n+1 
        do 3008 j = i, n+1 
          aa(1) = s(j, 1) - s(i-1, 1)
          b = aa(1)
          do 3010 k = 2, n 
            aa(k) = s(j, k) - s(i-1, k) - b
            b = b + aa(k)
 3010     continue

          call maxsubf(aa, m, fsum, mi1, mi2)
          if (fsum .gt. fmax) then
            fmax = fsum
            mind(1) = i-1
            mind(2) = j-1
            mind(3) = mi1
            mind(4) = mi2
          endif

 3008   continue
 3006 continue

      return
      end
