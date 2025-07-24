#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/species.h"
#include "fourdst/composition/composition.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"

#include "fourdst/config/config.h"

std::string EXAMPLE_FILENAME = std::string(getenv("MESON_SOURCE_ROOT")) + "/tests/config/example.yaml";

/**
 * @brief Test suite for the Composition class and related data structures.
 * @details This suite validates the core functionality of the composition library,
 * including the underlying atomic data, the creation and manipulation of
 * compositions, and the correctness of derived physical quantities.
 */
class compositionTest : public ::testing::Test {};

/**
 * @brief Tests the correctness of atomic mass data for select isotopes.
 * @details This test verifies that the static `fourdst::atomic::species` database is loaded
 * and that the `mass()` method for `Species` objects returns the expected values for a few
 * key isotopes (H-1, He-3, He-4).
 * @par What this test proves:
 * - The species data header is correctly included and parsed at compile time.
 * - The `species.at()` map lookup and the `mass()` accessor method work as expected for these specific cases.
 * @par What this test does not prove:
 * - The correctness of the mass for every isotope in the database. It is a spot check, not an exhaustive validation of the underlying data file.
 */
TEST_F(compositionTest, isotopeMasses) {
    EXPECT_NO_THROW(fourdst::atomic::species.at("H-1"));
    EXPECT_DOUBLE_EQ(fourdst::atomic::species.at("H-1").mass(), 1.007825031898);
    EXPECT_DOUBLE_EQ(fourdst::atomic::species.at("He-3").mass(), 3.0160293219700001);
    EXPECT_DOUBLE_EQ(fourdst::atomic::species.at("He-4").mass(),4.0026032541300003);
}

/**
 * @brief Tests the correctness of half-life data for select isotopes.
 * @details This test checks the `halfLife()` method for isotopes with different stability
 * characteristics: a stable isotope (H-1, infinite half-life), a radioactive isotope with a
 * finite half-life (F-18), and an unbound isotope (B-20, zero half-life).
 * @par What this test proves:
 * - The half-life data from the NUBASE data file is correctly parsed and stored.
 * - The `halfLife()` accessor correctly handles and returns values for stable (infinity), unstable (finite), and unbound (zero) cases.
 * @par What this test does not prove:
 * - The correctness of the half-life for all isotopes in the database.
 */
TEST_F(compositionTest, isotopeHalfLives) {
    EXPECT_DOUBLE_EQ(fourdst::atomic::H_1.halfLife(), std::numeric_limits<double>::infinity());
    EXPECT_DOUBLE_EQ(fourdst::atomic::F_18.halfLife(), 6584.04);
    EXPECT_DOUBLE_EQ(fourdst::atomic::B_20.halfLife(), 0.0);
}

/**
 * @brief Tests the numeric conversion of spin-parity strings.
 * @details This test validates the `spin()` method, which relies on the `convert_jpi_to_double`
 * utility. It covers a wide range of cases including half-integer spins (H-1), integer spins (Li-10),
 * zero spin (He-4), and cases where the spin is not known and should result in NaN (Bh-270).
 * @par What this test proves:
 * - The spin-parity string parsing logic correctly handles common formats (e.g., "1/2+", "5", "0+").
 * - The system correctly identifies and represents unknown or unmeasured spin values as `NaN`.
 * @par What this test does not prove:
 * - The correct parsing of every esoteric or malformed spin-parity string that might exist in data files. It focuses on the most common and expected formats.
 */
TEST_F(compositionTest, isotopeSpin) {
    using namespace fourdst::atomic;
    EXPECT_DOUBLE_EQ(H_1.spin(), 0.5);
    EXPECT_DOUBLE_EQ(He_4.spin(), 0.0);
    EXPECT_DOUBLE_EQ(Pm_164.spin(), 0.0);
    EXPECT_DOUBLE_EQ(Tb_164.spin(), 5.0);
    EXPECT_DOUBLE_EQ(Ta_163.spin(), 0.5);
    EXPECT_DOUBLE_EQ(Hf_165.spin(), 2.5);
    EXPECT_DOUBLE_EQ(Ta_165.spin(), 0.5);
    EXPECT_DOUBLE_EQ(Li_10.spin(), 1.0);
    EXPECT_DOUBLE_EQ(He_9.spin(), 0.5);
    EXPECT_DOUBLE_EQ(F_18.spin(), 0.0);
    EXPECT_DOUBLE_EQ(B_20.spin(), 1.0);
    EXPECT_TRUE(std::isnan(Bh_270.spin()));
}

/**
 * @brief Tests the default constructor of the Composition class.
 * @details This is a basic sanity check to ensure that a `Composition` object can be
 * instantiated without any arguments and does not throw an exception.
 * @par What this test proves:
 * - The default constructor is accessible and does not fail on basic initialization.
 * @par What this test does not prove:
 * - The state of the constructed object or the correctness of any of its methods.
 */
TEST_F(compositionTest, constructor) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    EXPECT_NO_THROW(fourdst::composition::Composition comp);
}

/**
 * @brief Tests the registration of valid and invalid symbols.
 * @details This test verifies that `registerSymbol` correctly adds valid isotope symbols
 * to the composition and throws an `InvalidSymbolError` for symbols that do not exist
 * in the atomic species database. It also checks that `getRegisteredSymbols` reflects the correct state.
 * @par What this test proves:
 * - The validation logic within `registerSymbol` correctly distinguishes between known and unknown species.
 * - The internal set of registered symbols is correctly maintained.
 * @par What this test does not prove:
 * - The handling of mode conflicts (e.g., trying to register a symbol in number fraction mode when the composition is already in mass fraction mode).
 */
TEST_F(compositionTest, registerSymbol) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    EXPECT_NO_THROW(comp.registerSymbol("H-1"));
    EXPECT_NO_THROW(comp.registerSymbol("He-4"));
    EXPECT_THROW(comp.registerSymbol("H-19"), fourdst::composition::exceptions::InvalidSymbolError);
    EXPECT_THROW(comp.registerSymbol("He-21"), fourdst::composition::exceptions::InvalidSymbolError);

    std::set<std::string> registeredSymbols = comp.getRegisteredSymbols();
    EXPECT_TRUE(registeredSymbols.find("H-1") != registeredSymbols.end());
    EXPECT_TRUE(registeredSymbols.find("He-4") != registeredSymbols.end());
    EXPECT_TRUE(registeredSymbols.find("H-19") == registeredSymbols.end());
    EXPECT_TRUE(registeredSymbols.find("He-21") == registeredSymbols.end());
}

/**
 * @brief Tests the core workflow of setting and getting mass fractions.
 * @details This test covers setting mass fractions for single and multiple symbols,
 * the requirement to `finalize()` before getting data, and the behavior when finalization
 * fails (e.g., due to non-normalized fractions).
 * @par What this test proves:
 * - `setMassFraction` correctly updates the internal state.
 * - `finalize` correctly validates the composition (sum of fractions must be ~1.0).
 * - `getComposition` and other getter methods correctly throw `CompositionNotFinalizedError` if called before `finalize`.
 * - An attempt to set a fraction for an unregistered symbol throws `UnregisteredSymbolError`.
 * @par What this test does not prove:
 * - The correctness of number fraction mode, which is tested separately.
 */
TEST_F(compositionTest, setGetComposition) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");

    EXPECT_DOUBLE_EQ(comp.setMassFraction("H-1", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setMassFraction("He-4", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setMassFraction("H-1", 0.6), 0.5);
    EXPECT_DOUBLE_EQ(comp.setMassFraction("He-4", 0.4), 0.5);

    EXPECT_NO_THROW(comp.finalize());
    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.6);

    EXPECT_THROW(comp.setMassFraction("He-3", 0.3), fourdst::composition::exceptions::UnregisteredSymbolError);

    const std::vector<std::string> symbols = {"H-1", "He-4"};

    EXPECT_NO_THROW(comp.setMassFraction(symbols, {0.5, 0.5}));
    EXPECT_THROW(auto r = comp.getComposition("H-1"), fourdst::composition::exceptions::CompositionNotFinalizedError);
    EXPECT_TRUE(comp.finalize());
    EXPECT_DOUBLE_EQ(comp.getComposition("H-1").first.mass_fraction(), 0.5);

    EXPECT_NO_THROW(comp.setMassFraction(symbols, {0.6, 0.6}));
    EXPECT_FALSE(comp.finalize());
    EXPECT_THROW(auto r = comp.getComposition("H-1"), fourdst::composition::exceptions::CompositionNotFinalizedError);
}

/**
 * @brief Tests the workflow of setting and getting number fractions.
 * @details This test mirrors `setGetComposition` but for number fraction mode. It verifies
 * that symbols can be registered in number fraction mode and that `setNumberFraction` and
 * `getNumberFraction` work as expected.
 * @par What this test proves:
 * - The composition can be correctly initialized and operated in number fraction mode.
 * - `setNumberFraction` and `getNumberFraction` function correctly.
 * - An attempt to set a fraction for an unregistered symbol throws the correct exception.
 * @par What this test does not prove:
 * - The correctness of conversions between mass and number fraction modes.
 */
TEST_F(compositionTest, setGetNumberFraction) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1", false);
    comp.registerSymbol("He-4", false);

    EXPECT_DOUBLE_EQ(comp.setNumberFraction("H-1", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setNumberFraction("He-4", 0.5), 0.0);
    EXPECT_DOUBLE_EQ(comp.setNumberFraction("H-1", 0.6), 0.5);
    EXPECT_DOUBLE_EQ(comp.setNumberFraction("He-4", 0.4), 0.5);

    EXPECT_NO_THROW(comp.finalize());
    EXPECT_DOUBLE_EQ(comp.getNumberFraction("H-1"), 0.6);

    EXPECT_THROW(comp.setNumberFraction("He-3", 0.3), fourdst::composition::exceptions::UnregisteredSymbolError);
}

/**
 * @brief Tests the creation of a normalized subset of a composition.
 * @details This test creates a composition, finalizes it, and then extracts a subset
 * containing only one of the original elements. It verifies that the `subset` method with
 * the "norm" option creates a new, valid composition where the single element's mass
 * fraction is normalized to 1.0.
 * @par What this test proves:
 * - The `subset` method can extract a subset of symbols.
 * - The "norm" method correctly renormalizes the fractions in the new subset to sum to 1.0.
 * @par What this test does not prove:
 * - The behavior of the `subset` method with the "none" option.
 */
TEST_F(compositionTest, subset) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6);
    comp.setMassFraction("He-4", 0.4);
    EXPECT_NO_THROW(comp.finalize());

    std::vector<std::string> symbols = {"H-1"};
    fourdst::composition::Composition subsetComp = comp.subset(symbols, "norm");
    EXPECT_TRUE(subsetComp.finalize());
    EXPECT_DOUBLE_EQ(subsetComp.getMassFraction("H-1"), 1.0);
}

/**
 * @brief Tests the auto-normalization feature of the `finalize` method.
 * @details This test sets mass fractions that do not sum to 1.0 and then calls
 * `finalize(true)`. It verifies that the composition is successfully finalized and that
 * the mass fractions are correctly scaled to sum to 1.0.
 * @par What this test proves:
 * - `finalize(true)` correctly calculates the sum of fractions and normalizes each entry.
 * - The resulting composition is valid and its normalized values can be retrieved.
 * @par What this test does not prove:
 * - The behavior of `finalize(false)`, which is tested separately.
 */
TEST_F(compositionTest, finalizeWithNormalization) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.3);
    comp.setMassFraction("He-4", 0.3);
    EXPECT_TRUE(comp.finalize(true));
    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.5);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.5);
}

/**
 * @brief Tests the default (non-normalizing) behavior of the `finalize` method.
 * @details This test sets mass fractions that already sum to 1.0 and calls `finalize(false)`.
 * It verifies that the composition is successfully finalized and the fractions remain unchanged.
 * @par What this test proves:
 * - `finalize(false)` or `finalize()` correctly validates a pre-normalized composition without altering its values.
 * @par What this test does not prove:
 * - That `finalize(false)` would fail for a non-normalized composition (this is implicitly tested in `setGetComposition`).
 */
TEST_F(compositionTest, finalizeWithoutNormalization) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.5);
    comp.setMassFraction("He-4", 0.5);
    EXPECT_TRUE(comp.finalize(false));
    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.5);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.5);
}

/**
 * @brief Tests the retrieval of global composition properties.
 * @details After creating and finalizing a composition, this test retrieves the
 * `CompositionEntry` and `GlobalComposition` data. It verifies that the mass fraction
 * in the entry and the calculated global properties (mean particle mass, specific number density)
 * are correct for the given input composition.
 * @par What this test proves:
 * - The `finalize` method correctly computes `meanParticleMass` and `specificNumberDensity`.
 * - The `getComposition` method returns a pair containing the correct entry-level and global data.
 * @par What this test does not prove:
 * - The correctness of these calculations for all possible compositions, particularly complex ones. It validates the mechanism for a simple binary mixture.
 */
TEST_F(compositionTest, getComposition) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6);
    comp.setMassFraction("He-4", 0.4);
    EXPECT_NO_THROW(comp.finalize());

    const auto compositionEntry = comp.getComposition("H-1");
    EXPECT_DOUBLE_EQ(compositionEntry.first.mass_fraction(), 0.6);
    EXPECT_DOUBLE_EQ(compositionEntry.second.meanParticleMass, 1.4382769310381101);
    EXPECT_DOUBLE_EQ(compositionEntry.second.specificNumberDensity, 1.0/1.4382769310381101);
}

/**
 * @brief Tests the ability to switch between mass and number fraction modes.
 * @details This test creates a composition in mass fraction mode, finalizes it, and then
 * switches to number fraction mode using `setCompositionMode(false)`. It then modifies the
 * composition using number fractions and verifies that it must be re-finalized before switching back.
 * @par What this test proves:
 * - `setCompositionMode` can be called on a finalized composition.
 * - After switching modes, the appropriate `set...Fraction` method can be used.
 * - Switching modes requires the composition to be finalized, and modifying it after the switch un-finalizes it again.
 * @par What this test does not prove:
 * - The numerical correctness of the fraction conversions that happen internally when the mode is switched.
 */
TEST_F(compositionTest, setCompositionMode) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
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

    EXPECT_THROW(comp.setCompositionMode(true), fourdst::composition::exceptions::CompositionNotFinalizedError);
    EXPECT_NO_THROW(comp.finalize());
    EXPECT_NO_THROW(comp.setCompositionMode(true));
}

/**
 * @brief Tests the `hasSymbol` utility method.
 * @details This test verifies that `hasSymbol` correctly reports the presence of registered
 * symbols and the absence of non-registered symbols.
 * @par What this test proves:
 * - The `hasSymbol` method accurately checks for the existence of a key in the internal composition map.
 * @par What this test does not prove:
 * - Anything about the state (e.g., mass fraction) of the symbol, only its presence.
 */
TEST_F(compositionTest, hasSymbol) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
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

/**
 * @brief Tests the mixing of two compositions.
 * @details This test creates two distinct compositions, finalizes them, and then mixes them
 * using both the `+` operator (50/50 mix) and the `mix` method with a specific fraction (25/75).
 * It verifies that the resulting mass fractions in the new compositions are correct.
 * @par What this test proves:
 * - The `mix` method and the `+` operator correctly perform linear interpolation of mass fractions.
 * - The resulting mixed composition is valid and its properties are correctly calculated.
 * @par What this test does not prove:
 * - The behavior when mixing compositions with non-overlapping sets of symbols.
 */
TEST_F(compositionTest, mix) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp1;
    comp1.registerSymbol("H-1");
    comp1.registerSymbol("He-4");
    comp1.setMassFraction("H-1", 0.6);
    comp1.setMassFraction("He-4", 0.4);
    EXPECT_NO_THROW(comp1.finalize());

    fourdst::composition::Composition comp2;
    comp2.registerSymbol("H-1");
    comp2.registerSymbol("He-4");
    comp2.setMassFraction("H-1", 0.4);
    comp2.setMassFraction("He-4", 0.6);
    EXPECT_NO_THROW(comp2.finalize());

    fourdst::composition::Composition mixedComp = comp1 + comp2;
    EXPECT_TRUE(mixedComp.finalize());
    EXPECT_DOUBLE_EQ(mixedComp.getMassFraction("H-1"), 0.5);
    EXPECT_DOUBLE_EQ(mixedComp.getMassFraction("He-4"), 0.5);

    fourdst::composition::Composition mixedComp2 = comp1.mix(comp2, 0.25);
    EXPECT_TRUE(mixedComp2.finalize());
    EXPECT_DOUBLE_EQ(mixedComp2.getMassFraction("H-1"), 0.45);
    EXPECT_DOUBLE_EQ(mixedComp2.getMassFraction("He-4"), 0.55);
}

/**
 * @brief Tests the calculation of molar abundance.
 * @details This test creates a simple composition and verifies that `getMolarAbundance`
 * returns the correct value, which is defined as (mass fraction / atomic mass).
 * @par What this test proves:
 * - The `getMolarAbundance` calculation is performed correctly.
 * @par What this test does not prove:
 * - The correctness of the underlying mass data, which is tested separately.
 */
TEST_F(compositionTest, molarAbundance) {
    fourdst::composition::Composition comp1;
    comp1.registerSymbol("H-1");
    comp1.registerSymbol("He-4");
    comp1.setMassFraction("H-1", 0.5);
    comp1.setMassFraction("He-4", 0.5);
    comp1.finalize();

    EXPECT_DOUBLE_EQ(comp1.getMolarAbundance("H-1"), 0.5/fourdst::atomic::H_1.mass());
    EXPECT_DOUBLE_EQ(comp1.getMolarAbundance("He-4"), 0.5/fourdst::atomic::He_4.mass());
}

/**
 * @brief Tests the registration and retrieval of species objects.
 * @details This test uses `registerSpecies` to add species directly (instead of by symbol)
 * and then uses `getRegisteredSpecies` to retrieve them. It verifies that the returned
 * set contains the correct `Species` objects.
 * @par What this test proves:
 * - The `registerSpecies` and `getRegisteredSpecies` methods work correctly with `Species` objects.
 * - The `std::set` of species correctly stores and orders the objects (based on the `<` operator overload).
 * @par What this test does not prove:
 * - The behavior of setting fractions via `Species` objects, which is covered in other tests.
 */
TEST_F(compositionTest, getRegisteredSpecies) {
    fourdst::composition::Composition comp;
    comp.registerSpecies({fourdst::atomic::Be_7, fourdst::atomic::H_1, fourdst::atomic::He_4}, true);
    auto registeredSpecies = comp.getRegisteredSpecies();
    EXPECT_TRUE(registeredSpecies.contains(fourdst::atomic::H_1));
    EXPECT_TRUE(registeredSpecies.contains(fourdst::atomic::He_4));
    EXPECT_FALSE(registeredSpecies.contains(fourdst::atomic::Li_6));
    auto it1 = registeredSpecies.begin();
    EXPECT_EQ(*it1, fourdst::atomic::H_1);
}

