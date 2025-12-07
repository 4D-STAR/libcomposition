#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <chrono>

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/composition.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "fourdst/composition/utils.h"
#include "fourdst/composition/decorators/composition_masked.h"
#include "fourdst/composition/utils/composition_hash.h"

#include "fourdst/config/config.h"


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
    fourdst::composition::Composition comp;
    EXPECT_NO_THROW(comp.registerSymbol("H-1"));
    EXPECT_NO_THROW(comp.registerSymbol("He-4"));
    EXPECT_THROW(comp.registerSymbol("H-19"), fourdst::composition::exceptions::UnknownSymbolError);
    EXPECT_THROW(comp.registerSymbol("He-21"), fourdst::composition::exceptions::UnknownSymbolError);

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
    // fourdst::config::Config::getInstance().loadConfig(EXAMPLE_FILENAME);
    fourdst::composition::Composition comp;
    EXPECT_NO_THROW(comp.registerSymbol("H-1"));
    EXPECT_NO_THROW(comp.registerSymbol("He-4"));

    EXPECT_NO_THROW(comp.setMolarAbundance("H-1", 0.6));
    EXPECT_NO_THROW(comp.setMolarAbundance("He-4", 0.4));

    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.27414655751871775);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.7258534424812823);

    EXPECT_THROW(comp.setMolarAbundance("He-3", 0.3), fourdst::composition::exceptions::UnregisteredSymbolError);

    EXPECT_NO_THROW(comp.registerSymbol("C-12"));
    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.27414655751871775);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.7258534424812823);

    EXPECT_NO_THROW(comp.setMolarAbundance("C-12", 0.1));

    EXPECT_DOUBLE_EQ(comp.getMassFraction("H-1"), 0.177551918933757);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("He-4"), 0.4701013674717613);
    EXPECT_DOUBLE_EQ(comp.getMassFraction("C-12"), 0.3523467135944818);

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
    comp.registerSpecies({fourdst::atomic::Be_7, fourdst::atomic::H_1, fourdst::atomic::He_4});
    auto registeredSpecies = comp.getRegisteredSpecies();
    EXPECT_TRUE(registeredSpecies.contains(fourdst::atomic::H_1));
    EXPECT_TRUE(registeredSpecies.contains(fourdst::atomic::He_4));
    EXPECT_FALSE(registeredSpecies.contains(fourdst::atomic::Li_6));
    const auto it1 = registeredSpecies.begin();
    EXPECT_EQ(*it1, fourdst::atomic::H_1);
}

/**
 * @brief Tests the az_to_species utility for correct species lookup and error handling.
 * @details This test checks that az_to_species returns the correct Species for valid (A,Z) pairs,
 * and the correct error codes for invalid atomic number or element number.
 * @par What this test proves:
 * - az_to_species correctly maps valid (A,Z) to Species.
 * - az_to_species returns appropriate error codes for invalid input.
 */
TEST_F(compositionTest, getSpeciesFromAZ) {
    EXPECT_EQ(fourdst::atomic::O_12, fourdst::atomic::az_to_species(12, 8));
    EXPECT_EQ(fourdst::atomic::SpeciesErrorType::SPECIES_SYMBOL_NOT_FOUND, fourdst::atomic::az_to_species(120, 38).error());
    EXPECT_EQ(fourdst::atomic::SpeciesErrorType::ELEMENT_SYMBOL_NOT_FOUND, fourdst::atomic::az_to_species(120, 500).error());
}


/**
 * @brief Tests mean atomic number and electron abundance calculations.
 * @details This test verifies that getElectronAbundance and getMeanAtomicNumber return correct values
 * based on the composition.
 * @par What this test proves:
 * - getElectronAbundance and getMeanAtomicNumber are calculated correctly.
 */
TEST_F(compositionTest, meanElectronAbundance) {
    using fourdst::atomic::species;
    using fourdst::composition::Composition;

    Composition comp;
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");

    comp.setMolarAbundance("H-1", 0.6);
    comp.setMolarAbundance("He-4", 0.4);

    const double expectedYe = 0.6 * species.at("H-1").z() + 0.4 * species.at("He-4").z();

    EXPECT_NEAR(comp.getElectronAbundance(), expectedYe, 1e-12);
}

TEST_F(compositionTest, buildFromMassFractionVector) {
    using fourdst::atomic::Species;
    using namespace fourdst::atomic;
    using fourdst::composition::Composition;

    const std::vector<Species> sVec = {H_1, Mg_24, He_4, C_12};
    const std::vector<double>  massFractions = {0.7, 0.01, 0.28, 0.01};

    const Composition comp = fourdst::composition::buildCompositionFromMassFractions(sVec, massFractions);

    EXPECT_DOUBLE_EQ(comp.getMassFraction(H_1), 0.7);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(He_4), 0.28);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(C_12), 0.01);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(Mg_24), 0.01);
}

TEST_F(compositionTest, buildFromMassFractionVectorString) {
    using fourdst::atomic::Species;
    using namespace fourdst::atomic;
    using fourdst::composition::Composition;

    const std::vector<std::string> sVec = {"H-1", "Mg-24", "He-4", "C-12"};
    const std::vector<double>  massFractions = {0.7, 0.01, 0.28, 0.01};

    const Composition comp = fourdst::composition::buildCompositionFromMassFractions(sVec, massFractions);

    EXPECT_DOUBLE_EQ(comp.getMassFraction(H_1), 0.7);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(He_4), 0.28);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(C_12), 0.01);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(Mg_24), 0.01);
}

TEST_F(compositionTest, buildFromMassFractionSet) {
    using fourdst::atomic::Species;
    using namespace fourdst::atomic;
    using fourdst::composition::Composition;

    const std::set<Species> sVec = {H_1, He_4, C_12};
    const std::vector<double>  massFractions = {0.7, 0.28, 0.02};

    const Composition comp = fourdst::composition::buildCompositionFromMassFractions(sVec, massFractions);

    EXPECT_DOUBLE_EQ(comp.getMassFraction(H_1), 0.7);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(He_4), 0.28);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(C_12), 0.02);
}

TEST_F(compositionTest, buildFromMassFractionUnorderedMap) {
    using fourdst::atomic::Species;
    using namespace fourdst::atomic;
    using fourdst::composition::Composition;

    const std::unordered_map<Species, double> sMap = {{H_1, 0.7}, {Mg_24, 0.01}, {He_4, 0.28}, {C_12, 0.01}};

    const Composition comp = fourdst::composition::buildCompositionFromMassFractions(sMap);

    EXPECT_DOUBLE_EQ(comp.getMassFraction(H_1), 0.7);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(He_4), 0.28);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(C_12), 0.01);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(Mg_24), 0.01);
}

TEST_F(compositionTest, buildFromMassFractionUnorderedMapString) {
    using fourdst::atomic::Species;
    using namespace fourdst::atomic;
    using fourdst::composition::Composition;

    const std::unordered_map<std::string, double> sMap = {{"H-1", 0.7}, {"Mg-24", 0.01}, {"He-4", 0.28}, {"C-12", 0.01}};

    const Composition comp = fourdst::composition::buildCompositionFromMassFractions(sMap);

    EXPECT_DOUBLE_EQ(comp.getMassFraction(H_1), 0.7);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(He_4), 0.28);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(C_12), 0.01);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(Mg_24), 0.01);
}

TEST_F(compositionTest, buildFromMassFractionOrderedMap) {
    using fourdst::atomic::Species;
    using namespace fourdst::atomic;
    using fourdst::composition::Composition;

    const std::map<Species, double> sMap = {{H_1, 0.7}, {Mg_24, 0.01}, {He_4, 0.28}, {C_12, 0.01}};

    const Composition comp = fourdst::composition::buildCompositionFromMassFractions(sMap);

    EXPECT_DOUBLE_EQ(comp.getMassFraction(H_1), 0.7);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(He_4), 0.28);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(C_12), 0.01);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(Mg_24), 0.01);
}

TEST_F(compositionTest, buildFromMassFractionOrderedMapString) {
    using fourdst::atomic::Species;
    using namespace fourdst::atomic;
    using fourdst::composition::Composition;

    const std::map<std::string, double> sMap = {{"H-1", 0.7}, {"Mg-24", 0.01}, {"He-4", 0.28}, {"C-12", 0.01}};

    const Composition comp = fourdst::composition::buildCompositionFromMassFractions(sMap);

    EXPECT_DOUBLE_EQ(comp.getMassFraction(H_1), 0.7);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(He_4), 0.28);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(C_12), 0.01);
    EXPECT_DOUBLE_EQ(comp.getMassFraction(Mg_24), 0.01);
}

TEST_F(compositionTest, decorators) {
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1"); comp.registerSymbol("He-4"); comp.registerSymbol("O-16");
    comp.setMolarAbundance("H-1", 0.6); comp.setMolarAbundance("He-4", 0.6);
    const fourdst::composition::MaskedComposition mComp(comp, {fourdst::atomic::H_1, fourdst::atomic::He_4});

    ASSERT_DOUBLE_EQ(mComp.getMolarAbundance(fourdst::atomic::H_1), 0.6);
    ASSERT_DOUBLE_EQ(mComp.getMolarAbundance("He-4"), 0.6);
    ASSERT_FALSE(mComp.contains("O-16"));

    comp.setMolarAbundance("H-1", 1.0);
    ASSERT_NE(mComp.getMolarAbundance(fourdst::atomic::H_1), 1.0);
}

TEST_F(compositionTest, orderInvariance) {
    fourdst::composition::Composition a;
    a.registerSymbol("He-4"); a.registerSymbol("H-1"); a.registerSymbol("O-16");
    a.setMolarAbundance("H-1", 0.6); a.setMolarAbundance("He-4", 0.6);

    fourdst::composition::Composition b;
    b.registerSymbol("O-16"); b.registerSymbol("H-1"); b.registerSymbol("He-4");
    b.setMolarAbundance("He-4", 0.6); b.setMolarAbundance("H-1", 0.6);

    const std::uint64_t ha = fourdst::composition::utils::CompositionHash::hash_exact(a);
    const std::uint64_t hb = fourdst::composition::utils::CompositionHash::hash_exact(b);
    ASSERT_EQ(ha, hb);
}

TEST_F(compositionTest, negativeZeroEqualsPositiveZero) {
    fourdst::composition::Composition a, b;
    a.registerSymbol("H-1"); b.registerSymbol("H-1");
    a.setMolarAbundance("H-1", 0.0);
    b.setMolarAbundance("H-1", -0.0);

    ASSERT_EQ(fourdst::composition::utils::CompositionHash::hash_exact(a),
              fourdst::composition::utils::CompositionHash::hash_exact(b));
}

TEST_F(compositionTest, cloneAndCopyStable) {
    std::vector<std::string> symbols = {"H-1", "He-4"};
    std::vector<double> abundances = {0.6, 0.4};
    fourdst::composition::Composition a(symbols, abundances);
    fourdst::composition::Composition b(a);

    const auto ha = fourdst::composition::utils::CompositionHash::hash_exact(a);
    const auto hb = fourdst::composition::utils::CompositionHash::hash_exact(b);
    ASSERT_EQ(ha, hb);

    std::unique_ptr<fourdst::composition::CompositionAbstract> cptr = a.clone();
    const auto hc = fourdst::composition::utils::CompositionHash::hash_exact(
        *static_cast<fourdst::composition::Composition*>(cptr.get()));
    ASSERT_EQ(ha, hc);
}

TEST_F(compositionTest, bothSidesRegisterSameZeroSpeciesEquality) {
    fourdst::composition::Composition a, b;
    a.registerSymbol("H-1"); b.registerSymbol("H-1");
    a.setMolarAbundance("H-1", 0.6); b.setMolarAbundance("H-1", 0.6);

    a.registerSymbol("He-4"); b.registerSymbol("He-4");
    ASSERT_EQ(fourdst::composition::utils::CompositionHash::hash_exact(a),
              fourdst::composition::utils::CompositionHash::hash_exact(b));
}

TEST_F(compositionTest, canonicalizeNaNIfAllowed) {
    fourdst::composition::Composition a, b;
    a.registerSymbol("H-1"); b.registerSymbol("H-1");
    double qnan1 = std::numeric_limits<double>::quiet_NaN();
    double qnan2 = std::bit_cast<double>(std::uint64_t{0x7ff80000'00000042ULL});
    a.setMolarAbundance("H-1", qnan1);
    b.setMolarAbundance("H-1", qnan2);
    ASSERT_EQ(fourdst::composition::utils::CompositionHash::hash_exact(a),
              fourdst::composition::utils::CompositionHash::hash_exact(b));
}

TEST_F(compositionTest, hash) {
    using namespace fourdst::atomic;
    fourdst::composition::Composition a, b;
    a.registerSymbol("H-1"); b.registerSymbol("H-1");
    a.setMolarAbundance("H-1", 0.6); b.setMolarAbundance("H-1", 0.6);
    EXPECT_NO_THROW((void)a.hash());
    EXPECT_EQ(a.hash(), b.hash());
    const std::unordered_map<Species, double> abundances ={
        {H_1, 0.702},
        {He_4, 0.06},
        {C_12, 0.001},
        {N_14, 0.0005},
        {O_16, 0.22},
    };
    fourdst::composition::Composition c(abundances), d(abundances);
    EXPECT_EQ(c.hash(), d.hash());
}

TEST_F(compositionTest, hashStaleing) {
    using namespace fourdst::atomic;
    const std::unordered_map<Species, double> abundances ={
        {H_1, 0.702},
        {He_4, 0.06},
        {C_12, 0.001},
        {N_14, 0.0005},
        {O_16, 0.22},
    };
    fourdst::composition::Composition a(abundances);

    const std::size_t hash1 = a.hash();
    a.setMolarAbundance("C-12", 0.002);
    const std::size_t hash2 = a.hash();
    EXPECT_NE(hash1, hash2);
}