// Entry point -- kept outside src/ so libskibidy.so has no main() symbol,
// allowing the test binary to provide its own gtest main.
#include "skibidy.h"

int main(int argc, const char** argv) {
  return bdm::skibidy::Simulate(argc, argv);
}
