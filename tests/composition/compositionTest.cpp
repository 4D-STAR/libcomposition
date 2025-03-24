#include <gtest/gtest.h>
#include <string>
#include <algorithm>

#include "atomicSpecies.h"

std::string EXAMPLE_FILENAME = std::string(getenv("MESON_SOURCE_ROOT")) + "/tests/composition/example.yaml";

/**
 * @brief Test suite for the composition class.
 */
class compositionTest : public ::testing::Test {};

/**
 * @brief Test the constructor of the composition class.
 */
TEST_F(compositionTest, constructor) {
    std::cout << "Testing the constructor of the composition class." << std::endl;
    std::cout << chemSpecies::species.at("H-1") << std::endl;
}

