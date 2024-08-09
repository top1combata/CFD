### Mesh concept

- face and cell indices are consecutive integers starts from 0

- each face has 2 neigbour cells, except for the boundary faces which have only 1 neighbour

- face's "owner" is the neighboring cell with the minimum index

- first element in every array in `m_face_neighbours` is always the owner cell index,  for the boundary face the second element of array is -1

- face vectors are directed outwards from the owner cell and have a length equal to its area

- in case of 2D "volumes" are areas and "areas" are lengths