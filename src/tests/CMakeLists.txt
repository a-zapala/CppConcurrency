add_executable(bottomlessBagTest bottomlessBagTest.cpp)
add_executable(sandArrangementTest sandArrangementTest.cpp)
add_executable(crystalSelectionTest crystalSelectionTest.cpp)


target_link_libraries( bottomlessBagTest pthread )

target_link_libraries( sandArrangementTest pthread )

target_link_libraries( crystalSelectionTest pthread )

