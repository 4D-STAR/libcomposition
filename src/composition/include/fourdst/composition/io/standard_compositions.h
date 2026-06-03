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
        AGS09,
        A09_Pryzbilla,
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
        {SolarCompositions::AGS09, "AGS09"},
        {SolarCompositions::A09_Pryzbilla, "A09_Pryzbilla"},
        {SolarCompositions::MB22_photospheric, "MB22_photospheric"},
        {SolarCompositions::AAG21_photospheric, "AAG21_photospheric"},
        {SolarCompositions::L09, "L09"}
    };

    inline std::unordered_map<IsotopicPercentages, std::string> IsotopicPercentages_to_string_map = {
        {IsotopicPercentages::L03, "L03_data"},
        {IsotopicPercentages::L09, "L09_data"}
    };

    /**
     * @class ChemicalFileParser
     * @brief An abstract base class for chemical file parsers.
     *
     * This class defines the interface for parsing fortran code files that contain
     * nuclide fractions. Derived classes must implement the `parse`
     * method to handle specific file formats.
     */
    class ChemicalFileParser {
    private:
        
    public:

        /**
         * @brief Parses a chemical file and returns the parsed data.
         *
         * This is a pure virtual function that must be implemented by derived
         * classes. It takes a filename as input and returns a `ParsedChemicalData`
         * struct containing the information extracted from the file.
         *
         * @param filename The path to the Chemical file to parse.
         * @return A `ParsedChemicalData` struct containing the parsed reaction data.
         *
         * @throws std::runtime_error If the file cannot be opened or a parsing
         * error occurs.
         *
         * @b Usage
         * @code
         * std::unique_ptr<ChemicalFileParser> parser = std::make_unique<SimpleReactionListFileParser>();
         * try {
         *     ParsedChemicalData data = parser->parse("my_reactions.txt");
         *     for (const auto& reaction_name : data.reactionPENames) {
         *         // ... process reaction name
        const mfem::GridFunction& grav_potential_at_inf(FEM& fem, const Args& args, const mfem::GridFunction& rho, bool pho_warm) {

}
         *     }
         * } catch (const std::runtime_error& e) {
         *     // ... handle error
         * }
         * @endcode
         */
        [[nodiscard]] CompositionData parse_compositon_data(const std::vector<char>& data,const std::string& scheme) const ;
        [[nodiscard]] IsotopicPercentage parse_isotopic_percentage(const std::vector<char>& data,const std::string& scheme) const ;
    };


}

namespace fourdst::composition {
    /**
     * @brief Function to retrieve a standard solar composition record indexed by their canonical names including
     *  - AG89
     *  - GN93
     *  - GS98
     *  - L03
     *  - AGS05
     *  - AGS08
     *  - A09_Pryzbilla
     *  - MB22_photospheric
     *  - AAG21_photospheric
     *  - L09
     *  Further, isotopic percentages can be selected as either
     *  - L03
     *  - L09
     *
     *  These data have been extracted from chem_def.f90 from MESA <version>
     *
     *  @note Composition names are case normalized; therefore, the inputs for metal fraction scheme and isotopic percentage scheme are case insensitive.
     *
     *  @param metal_fraction_scheme The name of the metal fraction scheme to use. Must be one of the following: AG89, GN93, GS98, L03, AGS05, AGS08, A09_Pryzbilla, MB22_photospheric, AAG21_photospheric, L09
     *  @param isotopic_percentage_scheme The name of the isotopic percentage scheme to use. Must be one of the following: L03, L09
     *  @param initial_z <poojan_documenent_here>
     *  @param initial_y <poojan document here>
     */
    [[nodiscard]] Composition get_composition_record(const std::string& metal_fraction_scheme,
                                                     const std::string& isotopic_percentage_scheme,
                                                     double initial_z,
                                                     double initial_y);

    /**
     * @brief Overload of the string based version of this function which accepts the enums Solar
     * @param metal_fraction_scheme Enum corresponding to the standard solar composition to select
     * @param isotopic_percentage_scheme Enum corresponding to the isotopic percentages prescription to select
     * @param initial_z <poojan_document_here>
     * @param initial_y <poojan_document_here>
     * @return
     */
    [[nodiscard]] Composition get_composition_record(io::SolarCompositions metal_fraction_scheme,
                                                     io::IsotopicPercentages isotopic_percentage_scheme,
                                                     double initial_z,
                                                     double initial_y);

}
