#include <iostream>
#include <memory>
#include <vector>

#include "../adventure.h"
#include "../utils.h"

void runAndVerify(Adventure &adventure, std::vector<Crystal> &crystals,
                  Crystal result) {
  Crystal crystal = adventure.selectBestCrystal(crystals);
  assert_msg(crystal == result, "Wrong crystal selection");
}

void testCase1(Adventure &adventure) {
  std::vector<Crystal> t1 = {Crystal(3), Crystal(2), Crystal(1)};
  Crystal r1 = Crystal(3);
  runAndVerify(adventure, t1, r1);
  std::vector<Crystal> t2 = {Crystal(5), Crystal(3), Crystal(2), Crystal(1),
                             Crystal(4)};
  Crystal r2 = Crystal(5);
  runAndVerify(adventure, t2, r2);
  std::vector<Crystal> t3 = {Crystal(7), Crystal(7), Crystal(7), Crystal(1),
                             Crystal(1), Crystal(4), Crystal(5)};
  Crystal r3 = Crystal(7);
  runAndVerify(adventure, t3, r3);
  std::vector<Crystal> t4 = {Crystal(1), Crystal(2), Crystal(3),
                             Crystal(4), Crystal(5), Crystal(6)};
  Crystal r4 = Crystal(6);
  runAndVerify(adventure, t4, r4);
  std::vector<Crystal> t5 = {Crystal(1), Crystal(2), Crystal(3), Crystal(9),
                             Crystal(4), Crystal(5), Crystal(6)};
  Crystal r5 = Crystal(9);
  runAndVerify(adventure, t5, r5);
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
      // });
    } else {
      std::vector<Crystal> t2(2575757);
      std::generate(t2.begin(), t2.end(), std::rand);
      Crystal r2 = *std::max_element(t2.begin(), t2.end());

      std::vector<Crystal> t3(1343434);
      std::generate(t3.begin(), t3.end(), std::rand);
      Crystal r3 = *std::max_element(t3.begin(), t3.end());

      // runAndPrintDuration( [&adventure, &t2, &r2]() {
      runAndVerify(*adventure, t2, r2);
      //});

      // runAndPrintDuration(    [&adventure, &t3, &r3]() {
      runAndVerify(*adventure, t3, r3);
      // });
    }
  }

  return 0;
}
