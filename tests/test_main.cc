// Explicit gtest main -- overrides main() from libskibidy.so
// (shared library symbols resolve after the binary's own symbols on Linux)
#include <gtest/gtest.h>

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
