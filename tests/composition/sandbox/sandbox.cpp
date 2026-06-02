#include "fourdst/composition/io/standard_compositions.h"
#include "fourdst/composition/io/StandardAbundancesBinary.h"

#include "fourdst/composition/composition.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <print>
#include <ranges>

int main(int argc, char** argv) {


    // @input: initial_z, initial_y, metal_fraction_scheme & isotopic_percentage_scheme
    // Options for metal_frac_scheme: ['AG89', 'GN93', 'GS98', 'L03', 'AGS05', 'AGSS09', 'A09_Przybilla', 'MB22_photospheric', 'AAG21_photospheric', 'L09']
    // Options for isotopic percentage scheme: [L03_data, L09_data]

    
    // CLI::App app("Loading Z fractions");

    // fourdst::config::Config<Options> config;
    // fourdst::config::register_as_cli(config, app);
    // app.parse(argc, argv)
    
    std::string metal_fraction_scheme, isotopic_percentage_scheme;
    double initial_z, initial_y;

    // the following four should be user input
    // initial_y can be optional
    initial_z = 0.02;
    initial_y = 0.24 + 2*initial_z;
    metal_fraction_scheme = "AG89";
    isotopic_percentage_scheme = "L03_data";

    fourdst::composition::io::ChemicalFileParser parser;
    fourdst::composition::Composition comp;
    comp = fourdst::composition::get_composition_record(metal_fraction_scheme, isotopic_percentage_scheme, initial_z, initial_y);
    std::cout << comp << std::endl;

}