The .gri file format is as follows:

nNode nElem Dim
for i = 1:nnode,
  x(i) y(i)
nBGroup
for i = 1:nBGroup
  nBFace(i) nf(i) Title(i)
  for j = 1:nBFace(i)
    B(i,j,1) .. B(i,j,nf(i))
nElem Order Basis
for i = 1:nElem
  E(i,1) .. E(i,nn)

where:
nNode     = number of nodes (vertices)
nElem     = number of elements
Dim       = dimension (2)
nBGroup   = number of boundary groups (i.e. BCs)
nBFace(i) = number of boundary edges in group i
nf(i)     = number of nodes per edge (2 in our case)
Title(i)  = title of boundary group i
NB(i,j,*) = list of boundary nodes for face j of group i
Order     = order of elements (1 in our case)
Basis     = basis of elements (triangles)
nn        = number of nodes per element (3)
E(i,*)    = list of nn nodes for element i
