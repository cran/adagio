subroutine maxsubf(x, n, s, i1, i2)
integer n, i1, i2, j1, j2
double precision x(n)
double precision s, s1, s2

s1 = 0; s2 = 0
i1 = 0; i2 = 0
j1 = 1; j2 = 1

do i = 1,n {
  if (s2 > -x(i)) {
    s2 = s2 + x(i)
    j2 = i
    if (s2 > s1) {
      s1 = s2
      i1 = j1; i2 = j2
    }
  } else {
    s2 = 0
    j1 = i+1; j2 = i+1
  }
}

return
end
