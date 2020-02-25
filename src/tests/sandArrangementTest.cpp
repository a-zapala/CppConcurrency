#include <vector>

#include "../adventure.h"
#include "../utils.h"

void runAndVerify(Adventure &adventure, std::vector<GrainOfSand> &grains,
                  std::vector<GrainOfSand> &result) {
  adventure.arrangeSand(grains);
  assert_msg(grains == result, "Wrong sand arrangement");
}

void testCase1(Adventure &adventure) {
  std::vector<GrainOfSand> t1 = {GrainOfSand(3), GrainOfSand(2),
                                 GrainOfSand(1)};
  std::vector<GrainOfSand> r1 = {GrainOfSand(1), GrainOfSand(2),
                                 GrainOfSand(3)};
  runAndVerify(adventure, t1, r1);
  std::vector<GrainOfSand> t2 = {GrainOfSand(5), GrainOfSand(3), GrainOfSand(2),
                                 GrainOfSand(1), GrainOfSand(4)};
  std::vector<GrainOfSand> r2 = {GrainOfSand(1), GrainOfSand(2), GrainOfSand(3),
                                 GrainOfSand(4), GrainOfSand(5)};
  runAndVerify(adventure, t2, r2);
  std::vector<GrainOfSand> t3 = {GrainOfSand(7), GrainOfSand(7), GrainOfSand(7),
                                 GrainOfSand(1), GrainOfSand(1), GrainOfSand(4),
                                 GrainOfSand(5)};
  std::vector<GrainOfSand> r3 = {GrainOfSand(1), GrainOfSand(1), GrainOfSand(4),
                                 GrainOfSand(5), GrainOfSand(7), GrainOfSand(7),
                                 GrainOfSand(7)};
  runAndVerify(adventure, t3, r3);
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
      //});
    } else {
      std::vector<GrainOfSand> t2(50000);
      std::generate(t2.begin(), t2.end(), std::rand);
      std::vector<GrainOfSand> r2 = t2;
      std::sort(r2.begin(), r2.end());

      std::vector<GrainOfSand> t3(11111);
      std::generate(t3.begin(), t3.end(), std::rand);
      std::vector<GrainOfSand> r3 = t3;
      std::sort(r3.begin(), r3.end());

      //    runAndPrintDuration(
      //      [&adventure, &t2, &r2]() {
      runAndVerify(*adventure, t2, r2);
      // });

      // runAndPrintDuration(    [&adventure, &t3, &r3]() {
      runAndVerify(*adventure, t3, r3);
      // });
    }
  }
  return 0;
}
