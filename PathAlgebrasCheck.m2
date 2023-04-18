TEST ///
-- test 0
-- needsPackage "PathAlgebras"
adj = matrix {{1,0},{1,1}}
G = paGraph({v,w},{e,f,g},adj)
R = QQ
A = R G
N = paMatrix {{e,f,g}}
I = {f*g - g^2}
M = A^1

Msyz1 = getSyzygies({M_0},I)  -- should stop with {}
assert(Msyz1 == {})
///

TEST ///
-- test 1
adj = matrix {{1,0},{1,1}}
G = paGraph({v,w},{e,f,g},adj)
R = QQ
A = R G
N = paMatrix {{e,f,g}}
I = {f*g - g^2}
M = A^1

Msyz2 = getSyzygies({M_0*v,M_0*w},I)  -- should have 2 syzygies
assert(Msyz2 / (p -> (entries p#0, p#1)) == {({w, 0}, {w, 0}), ({0, v}, {0, v})})
///

TEST ///
-- test 2
adj = matrix {{1,0},{1,1}}
G = paGraph({v,w},{e,f,g},adj)
R = QQ
A = R G
N = paMatrix {{e,f,g}}
I = {f*g - g^2}
M = A^1

Msyz3 = getSyzygies({M_0*g},I)
assert(Msyz3 / (p -> (entries p#0, p#1)) == {({v}, {v})})
///


TEST ///
-- test 3
adj = matrix {{1,0},{1,1}}
G = paGraph({v,w},{e,f,g},adj)
R = QQ
A = R G
N = paMatrix {{e,f,g}}
I = {f*g - g^2}
M = A^1

Msyz4 = getSyzygies({M_0*(e+f)},I)
assert(Msyz4 == {})
///


TEST ///
-- test 4
adj = matrix {{1,0},{1,1}}
G = paGraph({v,w},{e,f,g},adj)
R = QQ
A = R G
N = paMatrix {{e,f,g}}
I = {f*g - g^2}
M = A^1

Msyz5 = getSyzygies({M_0*e,M_0*f},I) -- should have 1 syzygy ???
assert(Msyz5 / (p -> (entries p#0, p#1)) == {({v, 0, 0}, {0, v}), ({0, w, 0}, {w, 0})})
///


TEST ///
-- test 5
adj = matrix {{1,0},{1,1}}
G = paGraph({v,w},{e,f,g},adj)
R = QQ
A = R G
N = paMatrix {{e,f,g}}
I = {f*g - g^2}
M = A^1

Msyz6 = getMinSyzygies({M_0*e,M_0*f},I)
assert(#Msyz6 == 2)
///

TEST ///
-- test 6
adj = matrix {{1,0},{1,1}}
G = paGraph({v,w},{e,f,g},adj)
R = QQ
A = R G
N = paMatrix {{e,f,g}}
I = {f*g - g^2}
M = A^1

Msyz7 = getSyzygies({M_0*e,M_0*e*f},I)
assert(Msyz7 / (p -> (entries p#0, p#1)) == {({w, 0},{w, 0}), ({-1*f, w + v},{-1*f, w + v})})
///

TEST ///
-- test 7
adj = matrix {{1,0},{1,1}}
G = paGraph({v,w},{e,f,g},adj)
R = QQ
A = R G
N = paMatrix {{e,f,g}}
I = {f*g - g^2}
M = A^1

Msyz8 = getSyzygies({M_0*e,2*M_0*e},I)
assert(Msyz8 / (p -> (entries p#0, p#1)) == {({w, 0},{0, 1/2*w}), ({w + v, -1/2*w + -1/2*v},{w + v, -1/2*w + -1/2*v})})
///

TEST ///
-- test 8
-- needsPackage "PathAlgebras"
adj = matrix {{3}}
G = paGraph({v},{x,y,z},adj)
R = QQ
A = R G
F0 = A^1
I = {x*y - y*x, x*z - z*x, y*z - z*y}
e0 = F0_0
f1 = e0*x
f2 = e0*y
f3 = e0*z
M = {f1,f2,f3}
M2 = getSyzygies(M,I)
assert( M2 / (p -> (entries p#0, p#1)) == {({-1*y, z, 0},{0, z, -1*y}), ({-1*x, 0, z},{z, 0, -1*x}), ({0, -1*x, y},{y, -1*x, 0})})
///


TEST ///
-- test 9
adj = matrix {{3}}
G = paGraph({v},{x,y,z},adj)
R = QQ
A = R G
F0 = A^1
I = {x*y - y*x, x*z - z*x, y*z - z*y}
e0 = F0_0
f1 = e0*x
f2 = e0*y
f3 = e0*z
M = {f1,f2,f3}
M2 = getSyzygies(M,I)
M3 = getSyzygies(M2,I)
assert(M3 / (p -> (entries p#0, p#1)) == {({x, -1*y, z},{x, -1*y, z})})
///

TEST ///
-- test 10
adj = matrix {{3}}
G = paGraph({v},{x,y,z},adj)
R = QQ
A = R G
F0 = A^1
I = {x*y - y*x, x*z - z*x, y*z - z*y}
e0 = F0_0
f1 = e0*x
f2 = e0*y
f3 = e0*z
M = {f1,f2,f3}
M2 = getSyzygies(M,I)
M3 = getSyzygies(M2,I)
M4 = getSyzygies(M3,I)
assert(M4 == {})
///





