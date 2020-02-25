#ifndef SRC_ADVENTURE_H_
#define SRC_ADVENTURE_H_

#include <algorithm>
#include <vector>

#include "../third_party/threadpool/threadpool.h"

#include "./types.h"
#include "./utils.h"

class Adventure {
 public:
  virtual ~Adventure() = default;

  virtual uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag &bag) = 0;

  virtual void arrangeSand(std::vector<GrainOfSand> &grains) = 0;

  virtual Crystal selectBestCrystal(std::vector<Crystal> &crystals) = 0;

 protected:
  static void merge(std::vector<GrainOfSand> &grains, uint64_t begin,
                    uint64_t mid, uint64_t end) {
    if (mid < end) {
      std::vector<GrainOfSand> pom;

      uint64_t i = begin;
      uint64_t j = mid;

      while (i < mid && j < end) {
        if (grains[i] < grains[j]) {
          pom.push_back(grains[i++]);
        } else {
          pom.push_back(grains[j++]);
        }
      }

      while (i < mid) {
        pom.push_back(grains[i++]);
      }

      while (j < end) {
        pom.push_back(grains[j++]);
      }

      for (uint64_t k = 0; k < pom.size(); k++) {
        grains[begin + k] = pom[k];
      }
    }
  }

  static void arrangeSandSequential(std::vector<GrainOfSand> &grains,
                                    uint64_t begin, uint16_t end) {
    uint64_t end_1 = begin + (end - begin) / 2;
    uint64_t end_2 = end;
    if (end - begin > 1) {
      arrangeSandSequential(grains, begin, end_1);
      arrangeSandSequential(grains, end_1, end_2);
      merge(grains, begin, end_1, end_2);
    }
  }

  static Crystal selectBestSequential(std::vector<Crystal> &crystals,
                                      uint64_t begin, uint64_t end) {
    uint64_t best_crystal = 0;

    for (uint64_t i = begin; i < end; i++) {
      if (crystals[best_crystal] < crystals[i]) best_crystal = i;
    }
    return crystals[best_crystal];
  }

  static std::vector<std::pair<uint64_t, uint64_t>> dividesInEqualPart(
      uint64_t size, uint64_t workers) {
    std::vector<std::pair<uint64_t, uint64_t>> result;

    uint64_t portion_size = size / workers;
    uint64_t current_begin = 0;

    for (uint64_t i = 0; i < workers; i++) {
      uint64_t current_end;
      if (i != workers - 1) {
        current_end = current_begin + portion_size;
      } else {
        current_end = size;
      }
      result.push_back(std::make_pair(current_begin, current_end));
      current_begin += portion_size;
    }

    return result;
  }
};

class LonesomeAdventure : public Adventure {
 public:
  LonesomeAdventure() {}

  virtual uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag &bag) {
    uint64_t bag_capacity = bag.getCapacity();
    uint64_t eggs_number = eggs.size();

    uint64_t **for_dynamic = new uint64_t *[eggs_number + 1];
    for (uint64_t i = 0; i <= eggs_number; ++i)
      for_dynamic[i] = new uint64_t[bag_capacity + 1];

    for (uint64_t j = 0; j <= bag_capacity; j++) {
      for_dynamic[0][j] = 0;
    }

    for (uint64_t i = 1; i <= eggs_number; i++) {
      for (uint64_t j = 0; j <= bag_capacity; j++) {
        if (eggs[i - 1].getSize() > j) {
          for_dynamic[i][j] = for_dynamic[i - 1][j];
        } else {
          for_dynamic[i][j] =
              std::max(for_dynamic[i - 1][j],
                       for_dynamic[i - 1][j - eggs[i - 1].getSize()] +
                           eggs[i - 1].getWeight());
        }
      }
    }

    uint64_t res = for_dynamic[eggs_number][bag_capacity];

    for (uint64_t i = 0; i <= eggs_number; ++i) delete[] for_dynamic[i];

    delete[] for_dynamic;
    return res;
  }

  virtual void arrangeSand(std::vector<GrainOfSand> &grains) {
    arrangeSandSequential(grains, 0, grains.size());
  }

  virtual Crystal selectBestCrystal(std::vector<Crystal> &crystals) {
    return selectBestSequential(crystals, 0, crystals.size());
  }
};

class TeamAdventure : public Adventure {
 private:
  class CyclicBarrier {
    std::condition_variable cv;
    std::mutex cv_m;
    uint64_t thread_number;
    unsigned index;
    uint64_t current_resistance[2];

   public:
    explicit CyclicBarrier(uint64_t resistance)
        : thread_number(resistance), index{0}, current_resistance{0, 0} {}

   public:
    void await() {
      std::unique_lock<std::mutex> lock(cv_m);
      unsigned current_index = index;
      current_resistance[current_index]++;

      if (current_resistance[current_index] < thread_number) {
        cv.wait(lock, [&] {
          return current_resistance[current_index] >= thread_number;
        });
      } else {
        index ^= 1;
        current_resistance[index] = 0;
        cv.notify_all();
      }
    }
  };

  static void packEggsThread(CyclicBarrier &barier, uint64_t begin_bag_capacity,
                             uint64_t end_bag_capacity, uint64_t **for_dynamic,
                             uint64_t eggs_number, std::vector<Egg> eggs,
                             const BottomlessBag &bag) {
    for (uint64_t i = 1; i <= eggs_number; i++) {
      for (uint64_t j = begin_bag_capacity; j < end_bag_capacity; j++) {
        if (eggs[i - 1].getSize() > j) {
          for_dynamic[i][j] = for_dynamic[i - 1][j];
        } else {
          for_dynamic[i][j] =
              std::max(for_dynamic[i - 1][j],
                       for_dynamic[i - 1][j - eggs[i - 1].getSize()] +
                           eggs[i - 1].getWeight());
        }
      }
      barier.await();
    }
  }

  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> getMergeTasks(
      std::vector<std::pair<uint64_t, uint64_t>> &equal_parts) {
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> result;
    for (uint64_t i = 0; i < equal_parts.size(); i += 2) {
      if (i != equal_parts.size() - 1) {
        result.push_back(std::make_tuple(equal_parts[i].first,
                                         equal_parts[i].second,
                                         equal_parts[i + 1].second));
      } else {
        result.push_back(std::make_tuple(equal_parts[i].first,
                                         equal_parts[i].second,
                                         equal_parts[i].second));
      }
    }
    return result;
  }

  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> joinMergeTask(
      const std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &merge_task) {
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> result;

    if (merge_task.size() == 1) return result;

    for (uint64_t i = 0; i < merge_task.size(); i += 2) {
      if (i != merge_task.size() - 1) {
        result.push_back(std::make_tuple(std::get<0>(merge_task[i]),
                                         std::get<0>(merge_task[i + 1]),
                                         std::get<2>(merge_task[i + 1])));
      } else {
        result.push_back(std::make_tuple(std::get<0>(merge_task[i]),
                                         std::get<2>(merge_task[i]),
                                         std::get<2>(merge_task[i])));
      }
    }
    return result;
  }

 public:
  explicit TeamAdventure(uint64_t numberOfShamansArg)
      : numberOfShamans(numberOfShamansArg),
        councilOfShamans(numberOfShamansArg) {}

  virtual void arrangeSand(std::vector<GrainOfSand> &grains) {
    auto equal_parts = dividesInEqualPart(grains.size(), numberOfShamans);
    std::future<void> *sorting_task = new std::future<void>[numberOfShamans];

    for (uint64_t i = 0; i < equal_parts.size(); i++) {
      sorting_task[i] = councilOfShamans.enqueue(
          TeamAdventure::arrangeSandSequential, std::ref(grains),
          equal_parts[i].first, equal_parts[i].second);
    }

    for (uint64_t i = 0; i < numberOfShamans; i++) {
      sorting_task[i].get();
    }
    delete[] sorting_task;

    auto merge_task_parts = getMergeTasks(equal_parts);
    while (merge_task_parts.size() > 0) {
      std::future<void> *merge_task =
          new std::future<void>[merge_task_parts.size()];
      for (uint64_t i = 0; i < merge_task_parts.size(); ++i) {
        merge_task[i] = councilOfShamans.enqueue(
            TeamAdventure::merge, std::ref(grains),
            std::get<0>(merge_task_parts[i]), std::get<1>(merge_task_parts[i]),
            std::get<2>(merge_task_parts[i]));
      }
      for (uint64_t i = 0; i < merge_task_parts.size(); i++) {
        merge_task[i].get();
      }
      delete[] merge_task;
      merge_task_parts = joinMergeTask(merge_task_parts);
    }
  }

  virtual uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag &bag) {
    uint64_t bag_capacity = bag.getCapacity();
    uint64_t eggs_number = eggs.size();

    uint64_t **for_dynamic = new uint64_t *[eggs_number + 1];

    for (uint64_t i = 0; i <= eggs_number; ++i)
      for_dynamic[i] = new uint64_t[bag_capacity + 1];

    for (uint64_t j = 0; j <= bag_capacity; j++) {
      for_dynamic[0][j] = 0;
    }

    auto equal_parts = dividesInEqualPart(bag_capacity + 1, numberOfShamans);

    std::future<void> *knapsack_tasks = new std::future<void>[numberOfShamans];
    CyclicBarrier barrier(numberOfShamans);

    for (uint64_t i = 0; i < equal_parts.size(); i++) {
      knapsack_tasks[i] = councilOfShamans.enqueue(
          TeamAdventure::packEggsThread, std::ref(barrier),
          equal_parts[i].first, equal_parts[i].second, for_dynamic, eggs_number,
          eggs, bag);
    }

    for (uint64_t i = 0; i < numberOfShamans; i++) {
      knapsack_tasks[i].get();
    }
    delete[] knapsack_tasks;
    uint64_t res = for_dynamic[eggs_number][bag_capacity];

    for (uint64_t i = 0; i <= eggs_number; ++i) delete[] for_dynamic[i];

    delete[] for_dynamic;
    return res;
  }

  virtual Crystal selectBestCrystal(std::vector<Crystal> &crystals) {
    auto equal_parts = dividesInEqualPart(crystals.size(), numberOfShamans);
    std::future<Crystal> *results = new std::future<Crystal>[numberOfShamans];

    for (uint64_t i = 0; i < equal_parts.size(); i++) {
      results[i] = councilOfShamans.enqueue(TeamAdventure::selectBestSequential,
                                            crystals, equal_parts[i].first,
                                            equal_parts[i].second);
    }

    Crystal best = results[0].get();
    for (uint64_t i = 1; i < numberOfShamans; i++) {
      Crystal best_in_partition = results[i].get();
      if (best < best_in_partition) best = best_in_partition;
    }
    delete[] results;
    return best;
  }

 private:
  uint64_t numberOfShamans;
  ThreadPool councilOfShamans;
};

#endif  // SRC_ADVENTURE_H_
