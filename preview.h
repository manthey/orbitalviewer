float CamXYZ[]={
   -7, -1,  1,    -7,  1,  1,    -7,  1, -1,    -7, -1, -1,
   -5, -1,  1,    -5,  1,  1,    -5,  1, -1,    -5, -1, -1,
   -3, -1,  1,    -3,  1,  1,    -3,  1, -1,    -3, -1, -1,
   -1, -1,  1,    -1,  1,  1,    -1,  1, -1,    -1, -1, -1,
    1, -1,  1,     1,  1,  1,     1,  1, -1,     1, -1, -1,
    3, -1,  1,     3,  1,  1,     3,  1, -1,     3, -1, -1,
    5, -1,  1,     5,  1,  1,     5,  1, -1,     5, -1, -1,
    7, -1,  1,     7,  1,  1,     7,  1, -1,     7, -1, -1,

    1, -7, -1,     1, -7,  1,    -1, -7,  1,    -1, -7, -1,
    1, -5, -1,     1, -5,  1,    -1, -5,  1,    -1, -5, -1,
    1, -3, -1,     1, -3,  1,    -1, -3,  1,    -1, -3, -1,
    1, -1, -1,     1, -1,  1,    -1, -1,  1,    -1, -1, -1,
    1,  1, -1,     1,  1,  1,    -1,  1,  1,    -1,  1, -1,
    1,  3, -1,     1,  3,  1,    -1,  3,  1,    -1,  3, -1,
    1,  5, -1,     1,  5,  1,    -1,  5,  1,    -1,  5, -1,
    1,  7, -1,     1,  7,  1,    -1,  7,  1,    -1,  7, -1,

   -1,  1, -7,     1,  1, -7,     1, -1, -7,    -1, -1, -7,
   -1,  1, -5,     1,  1, -5,     1, -1, -5,    -1, -1, -5,
   -1,  1, -3,     1,  1, -3,     1, -1, -3,    -1, -1, -3,
   -1,  1, -1,     1,  1, -1,     1, -1, -1,    -1, -1, -1,
   -1,  1,  1,     1,  1,  1,     1, -1,  1,    -1, -1,  1,
   -1,  1,  3,     1,  1,  3,     1, -1,  3,    -1, -1,  3,
   -1,  1,  5,     1,  1,  5,     1, -1,  5,    -1, -1,  5,
   -1,  1,  7,     1,  1,  7,     1, -1,  7,    -1, -1,  7};
short CamElem[26*15*3]={
   0, 1, 2, 3,  88, 1, 61, 1, 61,28, 88,28,  -1, 0, 0,

   0, 4, 5, 1,   1,28, 28,28, 28, 1,  1, 1,   0, 0, 1,
   1, 5, 6, 2,  28,61,  1,61,  1,88, 28,88,   0, 1, 0,
   2, 6, 7, 3,  28, 1,  1, 1,  1,28, 28,28,   0, 0,-1,
   3, 7, 4, 0,   1,88, 28,88, 28,61,  1,61,   0,-1, 0,
   4, 8, 9, 5,  31,88, 58,88, 58,61, 31,61,   0, 0, 1,
   5, 9,10, 6,  58, 1, 31, 1, 31,28, 58,28,   0, 1, 0,
   6,10,11, 7,  58,61, 31,61, 31,88, 58,88,   0, 0,-1,
   7,11, 8, 4,  31,28, 58,28, 58, 1, 31, 1,   0,-1, 0,
   8,12,13, 9,   1,28, 28,28, 28, 1,  1, 1,   0, 0, 1,
   9,13,14,10,  28,61,  1,61,  1,88, 28,88,   0, 1, 0,
  10,14,15,11,  28, 1,  1, 1,  1,28, 28,28,   0, 0,-1,
  11,15,12, 8,   1,88, 28,88, 28,61,  1,61,   0,-1, 0,

  16,20,21,17,   1,28, 28,28, 28, 1,  1, 1,   0, 0, 1,
  17,21,22,18, 118,31, 91,31, 91,58,118,58,   0, 1, 0,
  18,22,23,19,  28, 1,  1, 1,  1,28, 28,28,   0, 0,-1,
  19,23,20,16,  91,58,118,58,118,31, 91,31,   0,-1, 0,
  20,24,25,21, 121,58,148,58,148,31,121,31,   0, 0, 1,
  21,25,26,22,  58, 1, 31, 1, 31,28, 58,28,   0, 1, 0,
  22,26,27,23, 148,31,121,31,121,58,148,58,   0, 0,-1,
  23,27,24,20,  31,28, 58,28, 58, 1, 31, 1,   0,-1, 0,
  24,28,29,25,   1,28, 28,28, 28, 1,  1, 1,   0, 0, 1,
  25,29,30,26, 118,31, 91,31, 91,58,118,58,   0, 1, 0,
  26,30,31,27,  28, 1,  1, 1,  1,28, 28,28,   0, 0,-1,
  27,31,28,24,  91,58,118,58,118,31, 91,31,   0,-1, 0,

  29,28,31,30,  88, 1, 61, 1, 61,28, 88,28,   1, 0, 0};
