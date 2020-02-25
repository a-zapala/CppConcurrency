#! /usr/bin/env python
import collections
import os
import os.path
import subprocess

CLEANUP_PROJECT_CMD = (
    "git ls-files -o -i --exclude-standard --directory | xargs -r rm -r")

TESTS = {
    "bottomlessBagTest": [],
    "sandArrangementTest": [],
    "crystalSelectionTest": [],
    "bottomlessBagTest 1": [],
    "sandArrangementTest 1": [],
    "crystalSelectionTest 1": [],
}
PERFORMANCE_TESTS = [
    "bottomlessBagTest 1",
    "sandArrangementTest 1",
    "crystalSelectionTest 1",
]
SKIP_VALGRIND = [
    "bottomlessBagTest 1",
    "sandArrangementTest 1",
    "crystalSelectionTest 1",
]

PerformanceThreshold = collections.namedtuple(
    'PerformanceThreshold', ['name', 'threshold', 'less'])

PERFORMANCE_THRESHOLDS = []


TEST_RESULTS = {}
PERFORMANCE_RESULTS = {}
CHECKS = {}


def parse_test_output(test_name, output):
    res = {}
    for line in output.split('\n'):
        if line == "":
            continue
        tokens = line.split(';')
        if len(tokens) != 2:
            raise RuntimeError("Invalid output")
        metric = tokens[0]
        value = int(tokens[1])
        if metric not in TESTS[test_name]:
            raise RuntimeError("Unexpected metric (%s) in test %s" % (
                metric, test_name))
        if metric in res:
            raise RuntimeError("Duplicate metric (%s) in test %s" % (
                metric, test_name))
        res[metric] = value
    return res


def expect_no_output(cmd, error_msg):
    output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    if output == "":
        return True
    print error_msg
    print output
    return False


# Returns False if it found a memory leak, True otherwise.
def run_test(test, update_performance_stats, use_valgrind):
    test_path = os.path.join("src", "tests", test)
    output = None
    if use_valgrind:
        test_cmd = ("/usr/bin/valgrind --error-exitcode=123 "
                    "--leak-check=full %s") % test_path
    else:
        test_cmd = test_path
    try:
        output = subprocess.check_output("/usr/bin/timeout 180 " + test_cmd,
                                         shell=True)
        res = parse_test_output(test, output)
        if update_performance_stats:
            PERFORMANCE_RESULTS.update(res)
        TEST_RESULTS[test] = True
        print "== Test OK"
        return True
    except subprocess.CalledProcessError as ex:
        if use_valgrind and ex.returncode == 123:
            # It 's a memory leak. Don' t count it as a failure, only return
            # false
            return False
        print "== Error running test %s: %s" % (test, ex)
        print "== Test output:"
        print ex.output
        TEST_RESULTS[test] = False
        return True
    except Exception as ex:
        TEST_RESULTS[test] = False
        print "== Error running test %s: %s" % (test, ex)
        if output is not None:
            print "== Test output:"
            print output
        return True


def run_all_tests(update_performance_stats, use_valgrind):
    for test, _ in TESTS.iteritems():
        TEST_RESULTS[test] = TEST_RESULTS.setdefault(test, True)

        if not TEST_RESULTS[test]:
            print ("== Skipping test %s because it already failed in a "
                   "different config") % test
            continue

        if use_valgrind and test in SKIP_VALGRIND:
            print ("== Skipping test %s because it is not valgrind friendly"
                   % test)
            continue

        print "== Running test %s" % test
        res = run_test(test, update_performance_stats, use_valgrind)
        if use_valgrind and not res:
            print ("== Test %s discovered memory leaks, interrupting valgrind "
                   "checks") % test
            CHECKS["memory_leaks"] = False
            return


def was_test_fast_enough(test):
    for metric in TESTS[test]:
        for threshold in PERFORMANCE_THRESHOLDS:
            if threshold.name == metric:
                if threshold.less:
                    if PERFORMANCE_RESULTS[metric] > threshold.threshold:
                        print ("Performance problem: metric %s in test %s was "
                               "%s, which exceeded %s, performance -= 1") % (
                                   metric, test,
                                   PERFORMANCE_RESULTS[metric],
                                   threshold.threshold)
                        return False
                else:
                    if PERFORMANCE_RESULTS[metric] < threshold.threshold:
                        print ("Performance problem: metric %s in test %s "
                               "was %s, which didn't match %s, performance "
                               "-= 1") % (metric,
                                          test, PERFORMANCE_RESULTS[metric],
                                          threshold.threshold)
                        return False

    return True

def print_performance_results():
    print "=== Performance results"
    for test, _ in TESTS.iteritems():
        if test not in PERFORMANCE_TESTS:
            continue
        print "  test %s" % test
        for metric in TESTS[test]:
            if metric not in PERFORMANCE_RESULTS:
                continue
            for threshold in PERFORMANCE_THRESHOLDS:
                if threshold.name == metric:
                    print "    metric %s: %s (threshold %s)" % (
                        threshold.name, PERFORMANCE_RESULTS[metric],
                        threshold.threshold)


def summarize():
    correctness = 5
    code = 1
    performance = 2

    for test, _ in TESTS.iteritems():
        if not TEST_RESULTS[test]:
            correctness = correctness - 1
            print "Failed test %s, correctness -= 1" % test
    if not CHECKS["memory_leaks"]:
        correctness = correctness - 1
        print "Memory leaks, correctness -= 1"
    if correctness < 0:
        correctness = 0
    print "== Correctnes score: %s (max possible: 5)" % correctness

    for test, _ in TESTS.iteritems():
        if test not in PERFORMANCE_TESTS:
            continue
        if not TEST_RESULTS[test]:
            performance = performance - 1
            print ("Failed performance test %s, can't asses performance, "
                   "performance -= 1") % test
            continue
        if not was_test_fast_enough(test):
            performance = performance - 1
    if performance < 0:
        performance = 0
    print "== Performance score: %s (max possible: 2)" % performance

    if not CHECKS["lint"]:
        print "Linter warnings, code quality -= 1"
        code = code - 1
    if not CHECKS["formatting"]:
        print "Bad formatting, code quality -= 1"
        code = code - 1
#    if not CHECKS["clang-tidy"]:
#        print "clang-tidy discovered problems, code quality -= 1"
#        code = code - 1
    if performance + correctness < 4:
        print "Too few correctness/performance points to consider code quality"
        code = 0
    if code < 0:
        code = 0
    print "== Code quality score: %s (max possible: 1)" % code

    print "== Total score before reviewer adjustment: %s (max possible 8)" % (
        correctness + performance + code)

    print "== Reviewer's notes:"
    print
    print
    print
    print "== Final score: "


PROJECT_DIR = os.path.realpath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
print "== Changing directory to %s" % PROJECT_DIR
os.chdir(PROJECT_DIR)

print "== Cleaning project \"%s\"" % CLEANUP_PROJECT_CMD
subprocess.check_call(CLEANUP_PROJECT_CMD, shell=True,
                      stderr=subprocess.STDOUT)

print "== Configuring for debug"
subprocess.check_call(["cmake", "-DCMAKE_BUILD_TYPE=Debug"],
                      stderr=subprocess.STDOUT)

print "== Compiling"
subprocess.check_call(["make", "-j4"], stderr=subprocess.STDOUT)

print "== Running lint"
CHECKS["lint"] = (
    subprocess.call("./scripts/lint.sh", stderr=subprocess.STDOUT) == 0)

print "== Running format"
CHECKS["formatting"] = expect_no_output("./scripts/format.sh",
                                        "Code not properly formatted:")

print "== Running clang-tidy commented out"
#CHECKS["clang-tidy"] = expect_no_output("./scripts/tidy.sh",
#                                        "clang-tidy reported problems:")


print "========="
print "== Running all tests in Debug mode"
print "========="
run_all_tests(update_performance_stats=False, use_valgrind=False)

print "========="
print "== Running all tests in Debug mode with valgrind"
print "========="
CHECKS["memory_leaks"] = True
run_all_tests(update_performance_stats=False, use_valgrind=True)


print "== Cleaning project \"%s\"" % CLEANUP_PROJECT_CMD
subprocess.check_call(CLEANUP_PROJECT_CMD, shell=True,
                      stderr=subprocess.STDOUT)

print "== Configuring for release"
subprocess.check_call(["cmake", "-DCMAKE_BUILD_TYPE=Release"],
                      stderr=subprocess.STDOUT)

print "== Compiling"
subprocess.check_call(["make", "-j4"], stderr=subprocess.STDOUT)

print "========="
print "== Running all tests in Release mode"
print "========="
run_all_tests(update_performance_stats=True, use_valgrind=False)

print "== All done"

print "========="
print "== Summary (note that the real score will be calculated with different test set)"
print "========="

summarize()
