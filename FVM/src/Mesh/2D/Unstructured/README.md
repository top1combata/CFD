### PolyMesh format

```
# this is the comment
# the enumeration of all entities starts from 0

# vertices in (x,y) format
v 0 0
v 0 0.5
v 0.5 0
v 0.5 0.5

# faces in (v1, v2), where v1,v2 - indices of vertices

f 0 1
f 1 2
f 2 3
f 3 0

# cells in (f1, f2, f3, ...), where f1,f2,f3,... - indices of neighbour faces

c 0 1 2 3

# boundary conditions format
# b f p/U fixedValue/fixedGradient value
# f - face index, p or U boundary field
# fixedValue or fixedGradient - Dirichlet or Neumann boundary condition
# value - scalar in case of p, or (x,y) vector in case of U

b 0 p fixedGradient 0
b 0 U fixedValue 0 0
b 1 p fixedGradient 0
b 1 U fixedValue 0.25 0
...
```
