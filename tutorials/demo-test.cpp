#include "Snap.h"

int main(int argc, char* []) {
  TIntV vec;
  vec.Reserve(10, 10);
  printf("%d\n", vec.Len());
  vec[2] = 10;
  printf("%d\n", vec[2].Val);
  return 0;
}

