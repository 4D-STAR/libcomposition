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
 * - The state of the constructed object or the correctness of its methods.
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
    EXPECT_TRUE(registeredSymbols.contains("H-1"));
    EXPECT_TRUE(registeredSymbols.contains("He-4"));
    EXPECT_TRUE(!registeredSymbols.contains("H-19"));
    EXPECT_TRUE(!registeredSymbols.contains("He-21"));
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

    EXPECT_NO_THROW(static_cast<void>(comp.finalize()));
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

    EXPECT_NO_THROW(static_cast<void>(comp.finalize()));
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
    EXPECT_NO_THROW(static_cast<void>(comp.finalize()));

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
    EXPECT_NO_THROW(static_cast<void>(comp.finalize()));

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
    EXPECT_NO_THROW(static_cast<void>(comp.finalize()));

    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.6);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.4);

    EXPECT_NO_THROW(comp.setCompositionMode(false));

    EXPECT_NO_THROW(comp.setNumberFraction("H-1", 0.9));
    EXPECT_NO_THROW(comp.setNumberFraction("He-4", 0.1));

    EXPECT_THROW(comp.setCompositionMode(true), fourdst::composition::exceptions::CompositionNotFinalizedError);
    EXPECT_NO_THROW(static_cast<void>(comp.finalize()));
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
    EXPECT_NO_THROW(static_cast<void>(comp.finalize()));

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
    EXPECT_NO_THROW(static_cast<void>(comp1.finalize()));

    fourdst::composition::Composition comp2;
    comp2.registerSymbol("H-1");
    comp2.registerSymbol("He-4");
    comp2.setMassFraction("H-1", 0.4);
    comp2.setMassFraction("He-4", 0.6);
    EXPECT_NO_THROW(static_cast<void>(comp2.finalize()));

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
    const bool didFinalize = comp1.finalize();

    EXPECT_TRUE(didFinalize);
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

TEST_F(compositionTest, getSpeciesFromAZ) {
    EXPECT_EQ(fourdst::atomic::O_12, fourdst::atomic::az_to_species(12, 8));
    EXPECT_EQ(fourdst::atomic::SpeciesErrorType::SPECIES_SYMBOL_NOT_FOUND, fourdst::atomic::az_to_species(120, 38).error());
    EXPECT_EQ(fourdst::atomic::SpeciesErrorType::ELEMENT_SYMBOL_NOT_FOUND, fourdst::atomic::az_to_species(120, 500).error());
}

TEST_F(compositionTest, constructorWithSymbolsVectorAndSet) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    std::vector<std::string> vs = {"H-1", "He-4"};
    std::set<std::string> ss = {"H-1", "He-4"};

    fourdst::composition::Composition compVec(vs);
    EXPECT_TRUE(compVec.hasSymbol("H-1"));
    EXPECT_TRUE(compVec.hasSymbol("He-4"));

    fourdst::composition::Composition compSet(ss);
    EXPECT_TRUE(compSet.hasSymbol("H-1"));
    EXPECT_TRUE(compSet.hasSymbol("He-4"));
}

TEST_F(compositionTest, constructorWithSymbolsAndFractionsMassAndNumber) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;
    using fourdst::atomic::species;

    // Mass-fraction constructor
    std::vector<std::string> symM = {"H-1", "He-4"};
    std::vector<double> fracM = {0.6, 0.4};
    Composition compM(symM, fracM, true);
    EXPECT_NEAR(compM.getMassFraction("H-1"), 0.6, 1e-12);
    EXPECT_NEAR(compM.getMassFraction("He-4"), 0.4, 1e-12);
    // Mean particle mass and specific number density are reciprocals
    double sn = 0.6/species.at("H-1").mass() + 0.4/species.at("He-4").mass();
    double mp = 1.0/sn;
    EXPECT_NEAR(compM.getMeanParticleMass(), mp, 1e-12);

    // Number-fraction constructor
    std::vector<std::string> symN = {"H-1", "He-4"};
    std::vector<double> fracN = {0.9, 0.1};
    Composition compN(symN, fracN, false);
    EXPECT_NEAR(compN.getNumberFraction("H-1"), 0.9, 1e-12);
    EXPECT_NEAR(compN.getNumberFraction("He-4"), 0.1, 1e-12);
    double meanA = 0.9*species.at("H-1").mass() + 0.1*species.at("He-4").mass();
    EXPECT_NEAR(compN.getMeanParticleMass(), meanA, 1e-12);
    // Check converted mass fractions X_i = n_i * A_i / <A>
    double xH = 0.9*species.at("H-1").mass()/meanA;
    double xHe = 0.1*species.at("He-4").mass()/meanA;
    compN.setCompositionMode(true);
    EXPECT_NEAR(compN.getMassFraction("H-1"), xH, 1e-12);
    EXPECT_NEAR(compN.getMassFraction("He-4"), xHe, 1e-12);
}

TEST_F(compositionTest, registerSymbolVectorAndSingleSpeciesAndModeMismatch) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;
    Composition comp;
    comp.registerSymbol(std::vector<std::string>{"H-1", "He-4"});
    EXPECT_TRUE(comp.hasSymbol("H-1"));
    EXPECT_TRUE(comp.hasSymbol("He-4"));

    // Register by Species
    Composition comp2;
    comp2.registerSpecies(fourdst::atomic::H_1);
    comp2.registerSpecies(fourdst::atomic::He_4);
    EXPECT_TRUE(comp2.hasSymbol("H-1"));
    EXPECT_TRUE(comp2.hasSymbol("He-4"));

    // Mode mismatch should throw
    Composition comp3;
    comp3.registerSymbol("H-1", true); // mass mode
    EXPECT_THROW(comp3.registerSymbol("He-4", false), fourdst::composition::exceptions::CompositionModeError);
}

TEST_F(compositionTest, setMassFractionBySpeciesAndVector) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;
    using fourdst::atomic::H_1;
    using fourdst::atomic::He_4;

    Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");

    // Single species overload
    double old = comp.setMassFraction(H_1, 0.7);
    EXPECT_NEAR(old, 0.0, 1e-15);
    old = comp.setMassFraction(He_4, 0.3);
    EXPECT_NEAR(old, 0.0, 1e-15);

    // Vector overload
    std::vector<fourdst::atomic::Species> sp = {H_1, He_4};
    std::vector<double> xs = {0.6, 0.4};
    auto olds = comp.setMassFraction(sp, xs);
    ASSERT_EQ(olds.size(), 2u);
    EXPECT_NEAR(olds[0], 0.7, 1e-12);
    EXPECT_NEAR(olds[1], 0.3, 1e-12);

    EXPECT_TRUE(comp.finalize());
    EXPECT_NEAR(comp.getMassFraction("H-1"), 0.6, 1e-12);
    EXPECT_NEAR(comp.getMassFraction(He_4), 0.4, 1e-12);
}

TEST_F(compositionTest, setNumberFractionOverloads) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;
    using fourdst::atomic::H_1;
    using fourdst::atomic::He_4;

    Composition comp;
    comp.registerSymbol("H-1", false);
    comp.registerSymbol("He-4", false);

    // Single symbol
    double old = comp.setNumberFraction("H-1", 0.8);
    EXPECT_NEAR(old, 0.0, 1e-15);
    // Vector of symbols
    auto oldv = comp.setNumberFraction(std::vector<std::string>{"H-1", "He-4"}, std::vector<double>{0.75, 0.25});
    ASSERT_EQ(oldv.size(), 2u);
    EXPECT_NEAR(oldv[0], 0.8, 1e-12);
    EXPECT_NEAR(oldv[1], 0.0, 1e-12);

    // Species and vector<Species>
    old = comp.setNumberFraction(H_1, 0.7);
    EXPECT_NEAR(old, 0.75, 1e-12);
    auto oldsv = comp.setNumberFraction(std::vector<fourdst::atomic::Species>{H_1, He_4}, std::vector<double>{0.6, 0.4});
    ASSERT_EQ(oldsv.size(), 2u);
    EXPECT_NEAR(oldsv[0], 0.7, 1e-12);
    EXPECT_NEAR(oldsv[1], 0.25, 1e-12);

    EXPECT_TRUE(comp.finalize());
    EXPECT_NEAR(comp.getNumberFraction("H-1"), 0.6, 1e-12);
    EXPECT_NEAR(comp.getNumberFraction(He_4), 0.4, 1e-12);
}

TEST_F(compositionTest, mixErrorCases) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;

    Composition a; a.registerSymbol("H-1"); a.registerSymbol("He-4"); a.setMassFraction("H-1", 0.6); a.setMassFraction("He-4", 0.4);
    bool didFinalizeA = a.finalize();
    EXPECT_TRUE(didFinalizeA);
    Composition b; b.registerSymbol("H-1"); b.registerSymbol("He-4"); b.setMassFraction("H-1", 0.5); b.setMassFraction("He-4", 0.5);
    // Not finalized second comp
    EXPECT_THROW(static_cast<void>(a.mix(b, 0.5)), fourdst::composition::exceptions::CompositionNotFinalizedError);
    bool didFinalizeB = b.finalize();
    EXPECT_TRUE(didFinalizeB);
    // Invalid fraction
    EXPECT_THROW(static_cast<void>(a.mix(b, -0.1)), fourdst::composition::exceptions::InvalidCompositionError);
    EXPECT_THROW(static_cast<void>(a.mix(b, 1.1)), fourdst::composition::exceptions::InvalidCompositionError);
}

TEST_F(compositionTest, getMassAndNumberFractionMapsAndSpeciesOverloads) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;

    Composition comp; comp.registerSymbol("H-1"); comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6); comp.setMassFraction("He-4", 0.4);
    ASSERT_TRUE(comp.finalize());

    auto m = comp.getMassFraction();
    ASSERT_EQ(m.size(), 2u);
    EXPECT_NEAR(m.at("H-1"), 0.6, 1e-12);
    EXPECT_NEAR(m.at("He-4"), 0.4, 1e-12);
    EXPECT_NEAR(comp.getMassFraction(fourdst::atomic::H_1), 0.6, 1e-12);
    EXPECT_NEAR(comp.getMolarAbundance(fourdst::atomic::H_1), m.at("H-1")/fourdst::atomic::H_1.mass(), 1e-12);

    // Switch to number-fraction mode and verify number maps
    comp.setCompositionMode(false);
    // Must re-finalize after modifications (mode switch itself keeps values consistent but not finalized status changed? setCompositionMode requires to be finalized; here we just switched modes)
    // Set specific number fractions and finalize
    comp.setNumberFraction("H-1", 0.7);
    comp.setNumberFraction("He-4", 0.3);
    ASSERT_TRUE(comp.finalize());

    auto n = comp.getNumberFraction();
    ASSERT_EQ(n.size(), 2u);
    EXPECT_NEAR(n.at("H-1"), 0.7, 1e-12);
    EXPECT_NEAR(n.at("He-4"), 0.3, 1e-12);
    EXPECT_NEAR(comp.getNumberFraction(fourdst::atomic::He_4), 0.3, 1e-12);
}

TEST_F(compositionTest, meanAtomicNumberAndElectronAbundance) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::atomic::species;
    using fourdst::composition::Composition;

    Composition comp; comp.registerSymbol("H-1"); comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6); comp.setMassFraction("He-4", 0.4);
    ASSERT_TRUE(comp.finalize());

    // Compute expected Ye = sum(X_i * Z_i / A_i)
    constexpr double xH = 0.6, xHe = 0.4;
    const double aH = species.at("H-1").a();
    const double aHe = species.at("He-4").a();
    const double zH = species.at("H-1").z();
    const double zHe = species.at("He-4").z();
    const double expectedYe = xH*zH/aH + xHe*zHe/aHe;

    EXPECT_NEAR(comp.getElectronAbundance(), expectedYe, 1e-12);

    // <Z> = <A> * sum(X_i * Z_i / A_i)
    const double expectedZ = comp.getMeanParticleMass() * expectedYe;
    EXPECT_NEAR(comp.getMeanAtomicNumber(), expectedZ, 1e-12);
}

TEST_F(compositionTest, canonicalCompositionAndCaching) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;

    Composition comp; comp.registerSymbol("H-1"); comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6); comp.setMassFraction("He-4", 0.4);
    ASSERT_TRUE(comp.finalize());

    auto canon1 = comp.getCanonicalComposition();
    EXPECT_NEAR(canon1.X, 0.6, 1e-12);
    EXPECT_NEAR(canon1.Y, 0.4, 1e-12);
    EXPECT_NEAR(canon1.Z, 0.0, 1e-12);

    // Call again to exercise caching code path
    auto canon2 = comp.getCanonicalComposition();
    EXPECT_NEAR(canon2.X, 0.6, 1e-12);
    EXPECT_NEAR(canon2.Y, 0.4, 1e-12);
    EXPECT_NEAR(canon2.Z, 0.0, 1e-12);

    // Add a metal and re-check
    Composition comp2; comp2.registerSymbol("H-1"); comp2.registerSymbol("He-4"); comp2.registerSymbol("O-16");
    comp2.setMassFraction("H-1", 0.6); comp2.setMassFraction("He-4", 0.35); comp2.setMassFraction("O-16", 0.05);
    ASSERT_TRUE(comp2.finalize());
    auto canon3 = comp2.getCanonicalComposition(true);
    EXPECT_NEAR(canon3.X, 0.6, 1e-12);
    EXPECT_NEAR(canon3.Y, 0.35, 1e-12);
    EXPECT_NEAR(canon3.Z, 0.05, 1e-12);
}

TEST_F(compositionTest, vectorsAndIndexingAndSpeciesAtIndex) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;
    using fourdst::atomic::species;

    Composition comp; comp.registerSymbol("H-1"); comp.registerSymbol("He-4"); comp.registerSymbol("O-16");
    comp.setMassFraction("H-1", 0.5); comp.setMassFraction("He-4", 0.3); comp.setMassFraction("O-16", 0.2);
    ASSERT_TRUE(comp.finalize());

    // Mass fraction vector sorted by mass: H-1, He-4, O-16
    auto mv = comp.getMassFractionVector();
    ASSERT_EQ(mv.size(), 3u);
    EXPECT_NEAR(mv[0], 0.5, 1e-12);
    EXPECT_NEAR(mv[1], 0.3, 1e-12);
    EXPECT_NEAR(mv[2], 0.2, 1e-12);

    // Species indices ordering
    size_t iH = comp.getSpeciesIndex("H-1");
    size_t iHe = comp.getSpeciesIndex("He-4");
    size_t iO = comp.getSpeciesIndex("O-16");
    EXPECT_LT(iH, iHe);
    EXPECT_LT(iHe, iO);
    EXPECT_EQ(comp.getSpeciesIndex(fourdst::atomic::H_1), iH);
    EXPECT_EQ(comp.getSpeciesIndex(species.at("He-4")), iHe);

    // Species at index
    auto sAtHe = comp.getSpeciesAtIndex(iHe);
    EXPECT_EQ(std::string(sAtHe.name()), std::string("He-4"));

    // Number fraction vector after switching modes
    comp.setCompositionMode(false);
    // Tweak number fractions and finalize
    // Compute expected number fractions from original mass fractions first
    double denom = 0.5/species.at("H-1").mass() + 0.3/species.at("He-4").mass() + 0.2/species.at("O-16").mass();
    double nH_exp = (0.5/species.at("H-1").mass())/denom;
    double nHe_exp = (0.3/species.at("He-4").mass())/denom;
    double nO_exp = (0.2/species.at("O-16").mass())/denom;

    auto nv0 = comp.getNumberFractionVector();
    ASSERT_EQ(nv0.size(), 3u);
    EXPECT_NEAR(nv0[iH], nH_exp, 1e-12);
    EXPECT_NEAR(nv0[iHe], nHe_exp, 1e-12);
    EXPECT_NEAR(nv0[iO], nO_exp, 1e-12);

    // Molar abundance vector X_i/A_i in mass mode; switch back to mass mode to verify
    comp.setCompositionMode(true);
    bool didFinalize = comp.finalize(true);
    EXPECT_TRUE(didFinalize);
    auto av = comp.getMolarAbundanceVector();
    ASSERT_EQ(av.size(), 3u);
    EXPECT_NEAR(av[iH], 0.5/species.at("H-1").mass(), 1e-12);
    EXPECT_NEAR(av[iHe], 0.3/species.at("He-4").mass(), 1e-12);
    EXPECT_NEAR(av[iO], 0.2/species.at("O-16").mass(), 1e-12);
}

TEST_F(compositionTest, containsAndPreFinalizationGuards) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1"); comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6); comp.setMassFraction("He-4", 0.4);
    // contains should throw before finalize
    EXPECT_THROW(static_cast<void>(comp.contains(fourdst::atomic::H_1)), fourdst::composition::exceptions::CompositionNotFinalizedError);

    ASSERT_TRUE(comp.finalize());
    EXPECT_TRUE(comp.contains(fourdst::atomic::H_1));
    EXPECT_FALSE(comp.contains(fourdst::atomic::Li_6));
}

TEST_F(compositionTest, subsetNoneMethodAndNormalizationFlow) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1"); comp.registerSymbol("He-4"); comp.registerSymbol("O-16");
    comp.setMassFraction("H-1", 0.5); comp.setMassFraction("He-4", 0.3); comp.setMassFraction("O-16", 0.2);
    ASSERT_TRUE(comp.finalize());

    fourdst::composition::Composition sub = comp.subset(std::vector<std::string>{"H-1", "He-4"}, "none");
    // Not normalized: finalize without normalization should fail
    EXPECT_FALSE(sub.finalize(false));
    // With normalization, it should succeed and scale to sum to 1
    EXPECT_TRUE(sub.finalize(true));
    double sum = sub.getMassFraction("H-1") + sub.getMassFraction("He-4");
    EXPECT_NEAR(sum, 1.0, 1e-12);
    EXPECT_NEAR(sub.getMassFraction("H-1"), 0.5/(0.5+0.3), 1e-12);
    EXPECT_NEAR(sub.getMassFraction("He-4"), 0.3/(0.5+0.3), 1e-12);
}

TEST_F(compositionTest, copyConstructorAndAssignmentIndependence) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;

    Composition comp; comp.registerSymbol("H-1"); comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6); comp.setMassFraction("He-4", 0.4);
    ASSERT_TRUE(comp.finalize());

    Composition copy(comp); // copy ctor
    EXPECT_NEAR(copy.getMassFraction("H-1"), 0.6, 1e-12);
    EXPECT_NEAR(copy.getMassFraction("He-4"), 0.4, 1e-12);

    Composition assigned; assigned = comp; // assignment
    EXPECT_NEAR(assigned.getMassFraction("H-1"), 0.6, 1e-12);
    EXPECT_NEAR(assigned.getMassFraction("He-4"), 0.4, 1e-12);

    // Modify original and ensure copies do not change
    comp.setMassFraction("H-1", 0.7); comp.setMassFraction("He-4", 0.3); ASSERT_TRUE(comp.finalize());
    EXPECT_NEAR(copy.getMassFraction("H-1"), 0.6, 1e-12);
    EXPECT_NEAR(assigned.getMassFraction("He-4"), 0.4, 1e-12);
}

TEST_F(compositionTest, getCompositionBySpeciesAndAllEntries) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    using fourdst::composition::Composition;

    Composition comp; comp.registerSymbol("H-1"); comp.registerSymbol("He-4");
    comp.setMassFraction("H-1", 0.6); comp.setMassFraction("He-4", 0.4);
    ASSERT_TRUE(comp.finalize());

    auto pairBySpec = comp.getComposition(fourdst::atomic::H_1);
    EXPECT_NEAR(pairBySpec.first.mass_fraction(), 0.6, 1e-12);
    EXPECT_NEAR(pairBySpec.second.meanParticleMass, comp.getMeanParticleMass(), 1e-15);

    auto all = comp.getComposition();
    ASSERT_EQ(all.first.size(), 2u);
    EXPECT_NEAR(all.first.at("H-1").mass_fraction(), 0.6, 1e-12);
    EXPECT_NEAR(all.first.at("He-4").mass_fraction(), 0.4, 1e-12);
    EXPECT_NEAR(all.second.meanParticleMass, comp.getMeanParticleMass(), 1e-15);
}

TEST_F(compositionTest, iterationBeginEndAndIndexOutOfRange) {
    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1"); comp.registerSymbol("He-4"); comp.registerSymbol("O-16");
    comp.setMassFraction("H-1", 0.5); comp.setMassFraction("He-4", 0.3); comp.setMassFraction("O-16", 0.2);
    ASSERT_TRUE(comp.finalize());

    // Iterate and count entries
    size_t count = 0;
    for (auto it = comp.begin(); it != comp.end(); ++it) {
        count++;
    }
    EXPECT_EQ(count, 3u);

    // Out-of-range access
    EXPECT_THROW(static_cast<void>(comp.getSpeciesAtIndex(100)), std::out_of_range);
}

TEST_F(compositionTest, abstractBase) {
    class UnrestrictedComposition : public fourdst::composition::Composition {
    private:
        fourdst::atomic::Species m_species;
        const Composition& m_composition;
    public:
        UnrestrictedComposition(const Composition& base, const fourdst::atomic::Species& species):
        Composition(base),
        m_species(species),
        m_composition(base)
        {}

        double getMolarAbundance(const fourdst::atomic::Species &species) const override {
            if (species == m_species) {
                return 1.0;
            }
            return m_composition.getMolarAbundance(species);
        }
    };

    fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1"); comp.registerSymbol("He-4"); comp.registerSymbol("O-16");
    comp.setMassFraction("H-1", 0.5); comp.setMassFraction("He-4", 0.3); comp.setMassFraction("O-16", 0.2);
    ASSERT_TRUE(comp.finalize());

    const UnrestrictedComposition uComp(comp, fourdst::atomic::H_1);

    ASSERT_DOUBLE_EQ(uComp.getMolarAbundance(fourdst::atomic::H_1), 1.0);
    ASSERT_DOUBLE_EQ(uComp.getMassFraction("He-4"), comp.getMassFraction("He-4"));

}

