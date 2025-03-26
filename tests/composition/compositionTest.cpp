#include <gtest/gtest.h>
#include <stdexcept>
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

    EXPECT_DOUBLE_EQ(comp.setMassFraction("H-1", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setMassFraction("He-4", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setMassFraction("H-1", 0.6), 0.5);
    EXPECT_DOUBLE_EQ(comp.setMassFraction("He-4", 0.4), 0.5);

    EXPECT_NO_THROW(comp.finalize());
    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.6);

    EXPECT_THROW(comp.setMassFraction("He-3", 0.3), std::runtime_error);

    EXPECT_NO_THROW(comp.setMassFraction({"H-1", "He-4"}, {0.5, 0.5}));
    EXPECT_THROW(comp.getComposition("H-1"), std::runtime_error);
    EXPECT_TRUE(comp.finalize());
    EXPECT_DOUBLE_EQ(comp.getComposition("H-1").first.mass_fraction(), 0.5);

    EXPECT_NO_THROW(comp.setMassFraction({"H-1", "He-4"}, {0.6, 0.6}));
    EXPECT_FALSE(comp.finalize());
    EXPECT_THROW(comp.getComposition("H-1"), std::runtime_error);
}

TEST_F(compositionTest, setGetNumberFraction) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    comp.registerSymbol("H-1", false);
    comp.registerSymbol("He-4", false);

    EXPECT_DOUBLE_EQ(comp.setNumberFraction("H-1", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setNumberFraction("He-4", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setNumberFraction("H-1", 0.6), 0.5);
    EXPECT_DOUBLE_EQ(comp.setNumberFraction("He-4", 0.4), 0.5);

    EXPECT_NO_THROW(comp.finalize());
    EXPECT_DOUBLE_EQ(comp.getNumberFraction("H-1"), 0.6);

    EXPECT_THROW(comp.setNumberFraction("He-3", 0.3), std::runtime_error);
}

TEST_F(compositionTest, subset) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6);
    comp.setMassFraction("He-4", 0.4);
    EXPECT_NO_THROW(comp.finalize());

    std::vector<std::string> symbols = {"H-1"};
    composition::Composition subsetComp = comp.subset(symbols, "norm");
    EXPECT_TRUE(subsetComp.finalize());
    EXPECT_DOUBLE_EQ(subsetComp.getMassFraction("H-1"), 1.0);
}

TEST_F(compositionTest, finalizeWithNormalization) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.3);
    comp.setMassFraction("He-4", 0.3);
    EXPECT_TRUE(comp.finalize(true));
    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.5);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.5);
}

TEST_F(compositionTest, finalizeWithoutNormalization) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.5);
    comp.setMassFraction("He-4", 0.5);
    EXPECT_TRUE(comp.finalize(false));
    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.5);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.5);
}

TEST_F(compositionTest, getComposition) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6);
    comp.setMassFraction("He-4", 0.4);
    EXPECT_NO_THROW(comp.finalize());

    auto compositionEntry = comp.getComposition("H-1");
    EXPECT_DOUBLE_EQ(compositionEntry.first.mass_fraction(), 0.6);
    EXPECT_DOUBLE_EQ(compositionEntry.second.meanParticleMass, 1.4382769310381101);
    EXPECT_DOUBLE_EQ(compositionEntry.second.specificNumberDensity, 1.0/1.4382769310381101);
}

TEST_F(compositionTest, setCompositionMode) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6);
    comp.setMassFraction("He-4", 0.4);
    EXPECT_NO_THROW(comp.finalize());

    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.6);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.4);

    EXPECT_NO_THROW(comp.setCompositionMode(false));

    EXPECT_NO_THROW(comp.setNumberFraction("H-1", 0.9));
    EXPECT_NO_THROW(comp.setNumberFraction("He-4", 0.1));

    EXPECT_THROW(comp.setCompositionMode(true), std::runtime_error);
    EXPECT_NO_THROW(comp.finalize());
    EXPECT_NO_THROW(comp.setCompositionMode(true));
}

TEST_F(compositionTest, hasSymbol) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6);
    comp.setMassFraction("He-4", 0.4);
    EXPECT_NO_THROW(comp.finalize());

    EXPECT_TRUE(comp.hasSymbol("H-1"));
    EXPECT_TRUE(comp.hasSymbol("He-4"));
    EXPECT_FALSE(comp.hasSymbol("H-2"));
    EXPECT_FALSE(comp.hasSymbol("He-3"));
}

TEST_F(compositionTest, mix) {
    Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    composition::Composition comp1;
    comp1.registerSymbol("H-1");
    comp1.registerSymbol("He-4");
    comp1.setMassFraction("H-1", 0.6);
    comp1.setMassFraction("He-4", 0.4);
    EXPECT_NO_THROW(comp1.finalize());

    composition::Composition comp2;
    comp2.registerSymbol("H-1");
    comp2.registerSymbol("He-4");
    comp2.setMassFraction("H-1", 0.4);
    comp2.setMassFraction("He-4", 0.6);
    EXPECT_NO_THROW(comp2.finalize());

    composition::Composition mixedComp = comp1 + comp2;
    EXPECT_TRUE(mixedComp.finalize());
    EXPECT_DOUBLE_EQ(mixedComp.getMassFraction("H-1"), 0.5);
    EXPECT_DOUBLE_EQ(mixedComp.getMassFraction("He-4"), 0.5);

    composition::Composition mixedComp2 = comp1.mix(comp2, 0.25);
    EXPECT_TRUE(mixedComp2.finalize());
    EXPECT_DOUBLE_EQ(mixedComp2.getMassFraction("H-1"), 0.45);
    EXPECT_DOUBLE_EQ(mixedComp2.getMassFraction("He-4"), 0.55);
}