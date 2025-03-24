#include <gtest/gtest.h>
#include <string>
#include <algorithm>

#include "atomicSpecies.h"
#include "composition.h"
#include "config.h"

std::string EXAMPLE_FILENAME = std::string(getenv("MESON_SOURCE_ROOT")) + "/tests/config/example.yaml";

/**
 * @brief Test suite for the composition class.
 */
class compositionTest : public ::testing::Test {};

/**
 * @brief Test the constructor of the composition class.
 */
TEST_F(compositionTest, isotopeMasses) {
    EXPECT_NO_THROW(chemSpecies::species.at("H-1"));
    EXPECT_DOUBLE_EQ(chemSpecies::species.at("H-1").mass(), 1.007825031898);
    EXPECT_DOUBLE_EQ(chemSpecies::species.at("He-3").mass(), 3.0160293219700001);
    EXPECT_DOUBLE_EQ(chemSpecies::species.at("He-4").mass(),4.0026032541300003);
}

TEST_F(compositionTest, constructor) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    EXPECT_NO_THROW(composition::Composition comp);
}

TEST_F(compositionTest, registerSymbol) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    EXPECT_NO_THROW(comp.registerSymbol("H-1"));
    EXPECT_NO_THROW(comp.registerSymbol("He-4"));
    EXPECT_THROW(comp.registerSymbol("H-19"), std::runtime_error);
    EXPECT_THROW(comp.registerSymbol("He-21"), std::runtime_error);

    std::set<std::string> registeredSymbols = comp.getRegisteredSymbols();
    EXPECT_TRUE(registeredSymbols.find("H-1") != registeredSymbols.end());
    EXPECT_TRUE(registeredSymbols.find("He-4") != registeredSymbols.end());
    EXPECT_TRUE(registeredSymbols.find("H-19") == registeredSymbols.end());
    EXPECT_TRUE(registeredSymbols.find("He-21") == registeredSymbols.end());
}

TEST_F(compositionTest, setGetComposition) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");

    EXPECT_DOUBLE_EQ(comp.setComposition("H-1", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setComposition("He-4", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setComposition("H-1", 0.6), 0.5);
    EXPECT_DOUBLE_EQ(comp.setComposition("He-4", 0.4), 0.5);

    EXPECT_NO_THROW(comp.finalize());
    std::unordered_map<std::string, composition::CompositionEntry> compositions = comp.getCompositions();
    EXPECT_DOUBLE_EQ(compositions["H-1"].mass_fraction, 0.6);

    EXPECT_THROW(comp.setComposition("He-3", 0.3), std::runtime_error);

    EXPECT_NO_THROW(comp.setComposition({"H-1", "He-4"}, {0.5, 0.5}));
    EXPECT_THROW(comp.getComposition("H-1"), std::runtime_error);
    EXPECT_TRUE(comp.finalize());
    EXPECT_DOUBLE_EQ(comp.getComposition("H-1").mass_fraction, 0.5);

    EXPECT_NO_THROW(comp.setComposition({"H-1", "He-4"}, {0.6, 0.6}));
    EXPECT_FALSE(comp.finalize());
    EXPECT_THROW(comp.getComposition("H-1"), std::runtime_error);
}