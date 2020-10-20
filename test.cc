#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <sys/syscall.h>
#include <unistd.h>

int main(int argc, char* argv[]) {
  if(argc < 2) return - 1;
  std::cout << std::atoi(argv[1]) << std::endl;
  int b;
  char* bb;
  const int len(8);
  switch(std::atoi(argv[1])) {
  case 0:
    std::cout << ".text" << std::endl;
    for(int64_t i = 0; i < len; i ++)
      *(char*)((uint64_t)(void*)main + i) = 0x90;
    break;
  case 1:
    std::cout << "libc" << std::endl;
    for(int64_t i = 0; i < len; i ++)
      *(char*)((uint64_t)(void*)syscall + i) = 0x90;
    break;
  case 2:
    std::cout << "stack" << std::endl;
    for(int64_t i = 0; i < len; i ++)
      *(char*)((uint64_t)(void*)&b + i) = 0x90;
    break;
  case 3:
    std::cout << "heap" << std::endl;
#define HEAPL 60
    bb = (char*)malloc(HEAPL);
    for(int64_t i = 0; i < len; i ++)
      bb[HEAPL + i] = 0x90;
    free(bb);
    break;
  case 4:
    std::cout << "read random addr" << std::endl;
    bb = (char*)(uint64_t)((((uint64_t)arc4random()) << 32) + (uint64_t)arc4random());
    for(int64_t i = 0; i < len; i ++)
      std::cout << *(char*)(bb + i) << std::endl;
    break;
  case 5:
    std::cout << "read kernel .text" << std::endl;
    bb = (char*)0xffffffff90000000;
    for(int64_t i = 0; i < len; i ++)
      std::cout << *(char*)(bb + i) << std::endl;
    break;
  case 6:
    std::cout << "read .text" << std::endl;
    for(int64_t i = 0; i < len; i ++)
      std::cout << *(char*)((uint64_t)(void*)main + i) << std::endl;
    break;
  case 7:
    std::cout << "read libc" << std::endl;
    for(int64_t i = 0; i < len; i ++)
      std::cout << *(char*)((uint64_t)(void*)syscall + i) << std::endl;
    break;
  default:
    assert(0 && "no such method.");
  }
  std::cout << "skip through..." << std::endl;
  return 0;
}

