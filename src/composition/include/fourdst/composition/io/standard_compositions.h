#pragma once

#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"
#include "fourdst/composition/composition.h"

#include "quill/Logger.h"

#include <string>
#include <vector>

namespace fourdst::composition::io  {
    typedef std::vector<std::string> ParsedChemicalData;
    struct CompositionData {
        std::string comment_str;
        double he_abundance;
        bool requires_atomic_weight;
        std::vector<std::string> elements;
        std::vector<double> abundances;
    };
    struct IsotopicPercentage {
        std::string comment_str;
        std::vector<int> atomic_numbers;
        std::vector<std::string> elements;
        std::vector<int> mass_numbers;
        std::vector<double> percentages;
    };

    enum class SolarCompositions {
        AG89,
        GN93,
        GS98,
        L03,
        AGS05,
        AGSS09,
        A09_Przybilla,
        MB22_photospheric,
        AAG21_photospheric,
        L09
    };

    enum class IsotopicPercentages {
        L03,
        L09
    };

    inline std::unordered_map<SolarCompositions, std::string> SolarCompositions_to_string_map = {
        {SolarCompositions::AG89, "AG89"},
        {SolarCompositions::GN93, "GN93"},
        {SolarCompositions::GS98, "GS98"},
        {SolarCompositions::L03, "L03"},
        {SolarCompositions::AGS05, "AGS05"},
        {SolarCompositions::AGSS09, "AGSS09"},
        {SolarCompositions::A09_Przybilla, "A09_Przybilla"},
        {SolarCompositions::MB22_photospheric, "MB22_photospheric"},
        {SolarCompositions::AAG21_photospheric, "AAG21_photospheric"},
        {SolarCompositions::L09, "L09"}
    };

    inline std::unordered_map<IsotopicPercentages, std::string> IsotopicPercentages_to_string_map = {
        {IsotopicPercentages::L03, "L03_data"},
        {IsotopicPercentages::L09, "L09_data"}
    };

    /**
     * @brief Returns a read-only span over the embedded standard solar composition data.
     *
     * The data are stored as a compile-time binary array generated from the
     * embedded resource file at build time.  The returned span references static
     * storage directly — no copies are made.
     *
     * @return `std::span<const unsigned char>` Non-owning view over the raw tagged
     *         composition data bytes, suitable for passing to `ChemicalFileParser`.
     *
     * @par Examples
     * @code{.cpp}
     * auto raw = fourdst::composition::io::get_raw_standard_solar_composition_data();
     * std::vector<char> buf(raw.begin(), raw.end());
     * @endcode
     */
    std::span<const unsigned char> get_raw_standard_solar_composition_data();

    /**
     * @class ChemicalFileParser
     * @brief Parser for the tagged flat-text chemical composition data format.
     *
     * Reads data buffers that contain one or more named blocks delimited by
     * `BEGIN <scheme>` / `END <scheme>` sentinels and extracts either bulk metal
     * composition records or per-isotope percentage tables.
     */
    class ChemicalFileParser {
    private:

    public:

        /**
         * @brief Parses a named composition block from a tagged flat-text data buffer.
         *
         * Scans the buffer line by line for `BEGIN {scheme}`, then extracts five
         * fixed-position fields (comment, He abundance, atomic-weight flag, element
         * list, log10 abundance list) until `END {scheme}` is reached.  He abundance
         * and metal abundances are converted from log10 to linear scale via
         * `pow(10, x)` before storage.
         *
         * @param[in] data   Raw byte buffer (e.g., from
         *                   `get_raw_standard_solar_composition_data()`).
         * @param[in] scheme Block tag to extract (e.g., `"GS98"`, `"AGSS09"`).
         *
         * @return `CompositionData` Populated struct; default-initialized if the
         *         scheme is not found.
         *
         * @throws std::invalid_argument If a numeric field cannot be parsed by
         *         `std::stod`.
         * @throws std::out_of_range If a numeric field value is out of `double` range.
         *
         * @par Examples
         * @code{.cpp}
         * auto raw = fourdst::composition::io::get_raw_standard_solar_composition_data();
         * std::vector<char> buf(raw.begin(), raw.end());
         *
         * fourdst::composition::io::ChemicalFileParser parser;
         * auto comp = parser.parse_composition_data(buf, "GS98");
         *
         * for (size_t i = 0; i < comp.elements.size(); ++i) {
         *     std::println("{}: {:.6e}", comp.elements[i], comp.abundances[i]);
         * }
         * @endcode
         */
        [[nodiscard]] static CompositionData parse_composition_data(const std::vector<char>& data, const std::string& scheme);

        /**
         * @brief Parses a named isotopic-percentage block from a tagged flat-text data buffer.
         *
         * Scans the buffer for `BEGIN {scheme}`, then extracts five fixed-position
         * fields (comment, atomic numbers, element symbols, mass numbers, isotopic
         * percentages) until `END {scheme}` is reached.  Percentages are stored on
         * the 0-100 scale.
         *
         * @param[in] data   Raw byte buffer containing tagged isotopic percentage blocks.
         * @param[in] scheme Block tag to extract (e.g., `"L03_data"`, `"L09_data"`).
         *
         * @return `IsotopicPercentage` Populated struct; default-initialized if the
         *         scheme is not found.
         *
         * @throws std::invalid_argument If an integer or double field cannot be
         *         parsed by `std::stoi` / `std::stod`.
         * @throws std::out_of_range If any parsed value exceeds its target type range.
         *
         * @par Examples
         * @code{.cpp}
         * auto raw = fourdst::composition::io::get_raw_standard_solar_composition_data();
         * std::vector<char> buf(raw.begin(), raw.end());
         *
         * fourdst::composition::io::ChemicalFileParser parser;
         * auto iso = parser.parse_isotopic_percentage(buf, "L03_data");
         *
         * for (size_t i = 0; i < iso.elements.size(); ++i) {
         *     std::println("{}-{}: {:.4f}%",
         *                  iso.elements[i], iso.mass_numbers[i], iso.percentages[i]);
         * }
         * @endcode
         */
        [[nodiscard]] static IsotopicPercentage parse_isotopic_percentage(const std::vector<char>& data, const std::string& scheme);
    };


}

namespace fourdst::composition {
    /**
     * @brief Constructs a stellar `Composition` from a named solar metal-fraction
     *        scheme and isotopic-percentage table, scaled to the supplied bulk
     *        hydrogen and helium mass fractions.
     *
     * Available metal fraction schemes (extracted from MESA `chem_def.f90`):
     *  - `AG89`  (Anders & Grevesse 1989)
     *  - `GN93`  (Grevesse & Noels 1993)
     *  - `GS98`  (Grevesse & Sauval 1998)
     *  - `L03`   (Lodders 2003)
     *  - `AGS05` (Asplund, Grevesse & Sauval 2005)
     *  - `AGSS09` (Asplund et al. 2009)
     *  - `A09_Przybilla`
     *  - `MB22_photospheric`
     *  - `AAG21_photospheric`
     *  - `L09`   (Lodders 2009)
     *
     * Available isotopic percentage schemes:
     *  - `L03_data` (Lodders 2003)
     *  - `L09_data` (Lodders 2009)
     *
     * **Algorithm:**
     * 1. **Data loading** — The embedded binary `StandardMetalFractions` is copied
     *    into a `std::vector<char>` and parsed twice: once for `metal_fraction_scheme`
     *    and once for `isotopic_percentage_scheme`.
     * 2. **Species list** — The isotope table is iterated; any isotope whose element
     *    appears in the metals list or is `"H"` / `"He"` is looked up in the global
     *    `atomic::species` registry by `"<Element>-<A>"` and added to the list.
     * 3. **H and He mass fractions** — Four entries are prepended (H-1, H-2, He-3,
     *    He-4) using Anders & Grevesse (1989) solar He3/He4 ratio:
     *    - X(H-1) = clamp(1 - Z - Y, 0, 1)
     *    - X(H-2) = 0
     *    - X(He-3) = Y * xsol_He3 / (xsol_He3 + xsol_He4)
     *    - X(He-4) = Y * xsol_He4 / (xsol_He3 + xsol_He4)
     *    where xsol_He3 = 2.9291e-5 and xsol_He4 = 2.7521e-1.
     * 4. **Atomic-weight weighting** — When `CompositionData::requires_atomic_weight`
     *    is `true`, each metal's number-fraction abundance is multiplied by the
     *    atomic mass of its most-abundant isotope (determined from the isotopic table).
     * 5. **Normalisation** — Metal fractions are summed and normalised to unity.
     * 6. **Isotope distribution** — Per-isotope mass fractions are computed as:
     *    X_i = Z_total * f_E * (p_i * m_i) / sum_j(p_j * m_j)
     *    where f_E is the normalised metal fraction, p_i the isotopic percentage
     *    (0-100 scale), and m_i the isotope's atomic mass.
     * 7. **Renormalisation** — Metal mass fractions are rescaled so their sum
     *    equals Z_total exactly.
     * 8. **Assembly** — `buildCompositionFromMassFractions(species, massFracs)` builds
     *    the final `Composition` object.
     *
     * @param[in] metal_fraction_scheme      Block tag of the desired solar metal
     *            composition (e.g., `"GS98"`, `"AGSS09"`).  Case-sensitive; must
     *            match a `BEGIN`/`END` tag in the embedded data exactly.
     * @param[in] isotopic_percentage_scheme Block tag of the isotopic percentage
     *            table (e.g., `"L03_data"`, `"L09_data"`).
     * @param[in] initial_z                  Total metal mass fraction Z (0 <= Z < 1).
     * @param[in] initial_y                  Total helium mass fraction Y (0 <= Y < 1,
     *            with X + Y + Z <= 1 recommended).  X(H-1) is clamped to [0, 1]
     *            if the constraint is violated.
     *
     * @return `Composition` Fully populated composition object with per-isotope
     *         mass fractions normalised to `initial_z` and `initial_y`.
     *
     * @throws std::out_of_range If a species name derived from the isotopic table
     *         is absent from `atomic::species`, or if either scheme tag is not
     *         present in the embedded data.
     * @throws std::invalid_argument If numeric fields in the embedded data are
     *         malformed (propagated from `std::stod` / `std::stoi`).
     *
     * @par Examples
     * @code{.cpp}
     * // Grevesse & Sauval (1998) at Z = 0.02, Y = 0.28
     * fourdst::composition::Composition comp =
     *     fourdst::composition::get_composition_record("GS98", "L03_data", 0.02, 0.28);
     *
     * double x_h1 = comp.massFraction("H-1");
     * std::println("X(H-1) = {:.6f}", x_h1);  // approx 0.70
     * @endcode
     */
    [[nodiscard]] Composition get_composition_record(const std::string& metal_fraction_scheme,
                                                     const std::string& isotopic_percentage_scheme,
                                                     double initial_z,
                                                     double initial_y);

    /**
     * @brief Enum-based overload of `get_composition_record()`.
     *
     * Translates strongly-typed enum values to their canonical string
     * representations via `SolarCompositions_to_string_map` and
     * `IsotopicPercentages_to_string_map`, then delegates to the string-based
     * overload.  This overload is preferred as it prevents scheme name typos.
     *
     * @param[in] metal_fraction_scheme      Enum identifying the desired solar metal
     *            composition (e.g., `SolarCompositions::GS98`).
     * @param[in] isotopic_percentage_scheme Enum identifying the isotopic percentage
     *            table (e.g., `IsotopicPercentages::L03`).
     * @param[in] initial_z                  Total metal mass fraction Z (0 <= Z < 1).
     * @param[in] initial_y                  Total helium mass fraction Y (0 <= Y < 1).
     *
     * @return `Composition` Fully populated composition; see the string-based
     *         overload for the complete algorithm description.
     *
     * @throws std::out_of_range If the enum value is absent from its lookup map
     *         (should not occur with valid named enum members).
     *
     * @par Examples
     * @code{.cpp}
     * using namespace fourdst::composition;
     * using namespace fourdst::composition::io;
     *
     * // Asplund et al. (2009) at proto-solar Z and Y
     * Composition comp = get_composition_record(
     *     SolarCompositions::AGSS09,
     *     IsotopicPercentages::L09,
     *     0.0134, 0.2485
     * );
     *
     * double x_he4 = comp.massFraction("He-4");
     * std::println("X(He-4) = {:.6f}", x_he4);
     * @endcode
     */
    [[nodiscard]] Composition get_composition_record(io::SolarCompositions metal_fraction_scheme,
                                                     io::IsotopicPercentages isotopic_percentage_scheme,
                                                     double initial_z,
                                                     double initial_y);

}
