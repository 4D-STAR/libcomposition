#include "fourdst/composition/io/standard_compositions.h"
#include "fourdst/composition/composition.h"
#include "fourdst/atomic/species.h"

#include <string>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <ranges>
#include "CLI/CLI.hpp"

int main(int argc, char** argv) {


    // @input: initial_z, initial_y, metal_fraction_scheme & isotopic_percentage_scheme
    // Options for metal_frac_scheme: ['AG89', 'GN93', 'GS98', 'L03', 'AGS05', 'AGSS09', 'A09_Przybilla', 'MB22_photospheric', 'AAG21_photospheric', 'L09']
    // Options for isotopic percentage scheme: ['L03_data', 'L09_data']

    double initial_z;
    std::string metal_fraction_scheme;

    auto keys = fourdst::composition::io::SolarCompositions_to_string_map | std::views::values | std::ranges::to<std::vector>();
    
    CLI::App app("Example App To Load Solar Composition");
    app.add_option("-z,--initial_z", initial_z, "Initial Z")->required();
    app.add_option("-c,--metal-fraction-scheme", metal_fraction_scheme)->
        check(
            CLI::IsMember(
                keys,
                CLI::ignore_case)
        );

    CLI11_PARSE(app, argc, argv);
    
    std::string  isotopic_percentage_scheme;
    double initial_y;

    // the following four should be user input
    // initial_y can be optional
    initial_y = 0.24 + 2*initial_z;
    isotopic_percentage_scheme = "L03_data";

    fourdst::composition::Composition comp;
    comp = fourdst::composition::get_composition_record(metal_fraction_scheme, isotopic_percentage_scheme, initial_z, initial_y);
    std::cout << comp << std::endl;

}