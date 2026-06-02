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
    [[nodiscard]] Composition get_composition_record(const std::string& metal_fraction_scheme,
                                                                    const std::string& isotopic_percentage_scheme,
                                                                    double initial_z, double initial_y);

}