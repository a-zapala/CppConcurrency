#include <iostream>

#include "../adventure.h"
#include "../utils.h"

void correctnessTest(std::vector<Egg> eggs, BottomlessBag bag,
                     uint64_t expectedResults, Adventure &adventure) {
  uint64_t result = adventure.packEggs(eggs, bag);
  assert_eq_msg(result, expectedResults, "Unexpected packing result");
}

void testCase1(Adventure &adventure) {
  std::vector<Egg> eggs1{Egg(1, 1), Egg(2, 2), Egg(3, 3)};
  for (int i = 0; i < 10; ++i) {
    correctnessTest(eggs1, BottomlessBag(i), std::min(i, 6), adventure);
  }
}

void testCase2(Adventure &adventure) {
  std::vector<Egg> eggs2{Egg(5, 99999), Egg(1, 1), Egg(2, 2), Egg(3, 3),
                         Egg(1, 99999)};
  correctnessTest(eggs2, BottomlessBag(1), 99999, adventure);
  correctnessTest(eggs2, BottomlessBag(3), 99999 + 2, adventure);
  correctnessTest(eggs2, BottomlessBag(5), 99999 + 4, adventure);
  correctnessTest(eggs2, BottomlessBag(6), 2 * 99999, adventure);
}

void testCase3(Adventure &adventure) {
  std::vector<Egg> eggs;
  for (int i = 0; i < 33; ++i) {
    eggs.push_back(Egg(i, i * i + 7));
  }

  correctnessTest(eggs, BottomlessBag(100), 2969, adventure);
}

void testCase4(Adventure &adventure) {
  std::vector<Egg> eggs;
  for (int i = 0; i < 100; ++i) {
    eggs.push_back(Egg(i % 10, i * 3 + 2));
  }

  correctnessTest(eggs, BottomlessBag(10000), 15050, adventure);
}

// This test may not parallelize well. Why?
void testCase5(Adventure &adventure) {
  std::vector<Egg> eggs;
  for (int i = 0; i < 700; ++i) {
    eggs.push_back(Egg(i, i * 5 + 33));
  }

  correctnessTest(eggs, BottomlessBag(2000), 12079, adventure);
}

int main(int argc, char **argv) {
  for (std::shared_ptr<Adventure> adventure :
       std::vector<std::shared_ptr<Adventure> >{
           std::shared_ptr<Adventure>(new LonesomeAdventure{}),
           std::shared_ptr<Adventure>(new TeamAdventure(1)),
           std::shared_ptr<Adventure>(new TeamAdventure(2)),
           std::shared_ptr<Adventure>(new TeamAdventure(3)),
           std::shared_ptr<Adventure>(new TeamAdventure(4)),
           std::shared_ptr<Adventure>(new TeamAdventure(8))}) {
    if (argc == 1) {
      // runAndPrintDuration([&adventure]() {
      testCase1(*adventure);
      testCase2(*adventure);
      testCase3(*adventure);
      // });
    } else {
      // runAndPrintDuration([&adventure]() {
      testCase4(*adventure);
      //});

      // runAndPrintDuration([&adventure]() {
      testCase5(*adventure);
      //});
    }
  }
  return 0;
}
