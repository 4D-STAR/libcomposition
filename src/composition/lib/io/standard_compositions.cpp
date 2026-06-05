#include "fourdst/composition/io/standard_compositions.h"
#include "fourdst/composition/io/StandardMetalFractionsBinary.h"

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
#include <cctype>

namespace fourdst:: composition::io {
    namespace {
        inline void ltrim(std::string &s) {
            s.erase(
                s.begin(),
                std::ranges::find_if(s,
                                     [](const unsigned char ch) {
                                         return !std::isspace(ch);
                                     })
             );
        }

        inline void rtrim(std::string &s) {
            s.erase(
                std::find_if(
                    s.rbegin(),
                    s.rend(),
                    [](const unsigned char ch) {
                        return !std::isspace(ch);
                    }).base(),
                    s.end()
                );
        }

        inline void trim(std::string &s) {
            ltrim(s);
            rtrim(s);
        }


    }

    bool to_bool(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(),
                    [](unsigned char c){ return std::tolower(c); });

        return s == "true";
    }

    CompositionData ChemicalFileParser::parse_compositon_data(const std::vector<char>& data,const std::string& scheme) const {

        // get file and metal_fraction_scheme
        // Load the file
        // find the metal_fraction_scheme
        // return abundances

        // LOG_TRACE_L1(m_logger, "Parsing chemical abundance for: {}", scheme);

        bool debug = false;
        if (debug){
            std::println("Parsing chemical abundance for: {}", scheme);
        }

        std::istringstream stream(std::string(data.begin(), data.end()));

        // add error message if something goes wrong

        std::string line;
        int start_line = 0;
        int i = 0;

        CompositionData comp;

        while (std::getline(stream, line)) {

            // find where the end of the scheme block is
            auto end_pos = std::ranges::search(line,std::format("END {}", scheme));

            // exit if have reached the end of block
            if (!end_pos.empty()) {
                break;
            }

            if (start_line>0){
                const size_t colon_pos = line.find(':');
                line = line.substr(colon_pos+1);

                line.erase(std::remove_if(line.begin(), line.end(),
                    [](char c){ return c == '[' || c == ']'; }),
                    line.end());

                trim(line);
                std::string item;
                std::stringstream ss(line);
                double val;
                switch(i-start_line){
                    case 1:
                        comp.comment_str = line;
                        break;
                    case 2:
                        comp.he_abundance = std::pow(10.0,std::stod(line));
                        break;
                    case 3:
                        comp.requires_atomic_weight = to_bool(line);
                        break;
                    case 4:
                        while(std::getline(ss, item, ',')) {
                            comp.elements.push_back(item);
                        }
                        break;
                    case 5:
                        while(std::getline(ss, item, ',')) {
                            val = std::pow(10.0,std::stod(item));
                            comp.abundances.push_back(val);
                        }
                        break;
                }
            }

            // find where the start of the scheme block is
            auto start_pos = std::ranges::search(line, std::format("BEGIN {}", scheme));

            if (!start_pos.empty()) {
                start_line = i;
            }
            i+=1;
        }

        // if (start_pos==0):
        //     raise error ("Scheme {} not found", scheme)


        if (debug){
            std::println("he_abundance: {}", comp.he_abundance);
            std::println("requires_atomic_weight: {}", comp.requires_atomic_weight);
            std::println("elements: {}",comp.elements);
            std::println("abundances: {}", comp.abundances);
        }

        return comp;

    }

    IsotopicPercentage ChemicalFileParser::parse_isotopic_percentage(const std::vector<char>& data,const std::string& scheme) const {

        // get file and iso_scheme
        // Load the file
        // find the iso_scheme
        // get iso_comp data
        // IsotopicPercentage object

        bool debug = false;
        if (debug){
            std::println("Parsing Isotopic Percentage for: {}", scheme);
        }

        std::istringstream stream(std::string(data.begin(), data.end()));

        // add error message if something goes wrong
        ParsedChemicalData parsed;

        std::string line;
        int start_line = 0;
        int i = 0;

        IsotopicPercentage iso;

        while (std::getline(stream, line)) {

            // find where the end of the scheme block is
            auto end_pos = std::ranges::search(line,std::format("END {}", scheme));

            // exit if have reached the end of block
            if (!end_pos.empty()) {
                break;
            }

            if (start_line>0){
                const size_t colon_pos = line.find(':');
                line = line.substr(colon_pos+1);

                line.erase(std::remove_if(line.begin(), line.end(),
                    [](char c){ return c == '[' || c == ']'; }),
                    line.end());

                trim(line);
                std::string item;
                std::stringstream ss(line);
                parsed.push_back(line);
                switch(i-start_line){
                    case 1:
                        iso.comment_str = line;
                        break;
                    case 2:
                        while(std::getline(ss, item, ',')) {
                            iso.atomic_numbers.push_back(std::stoi(item));
                        }
                        break;
                    case 3:
                        while(std::getline(ss, item, ',')) {
                            iso.elements.push_back(item);
                        }
                        break;
                    case 4:
                        while(std::getline(ss, item, ',')) {
                            iso.mass_numbers.push_back(std::stoi(item));
                        }
                        break;
                    case 5:
                        while(std::getline(ss, item, ',')) {
                            iso.percentages.push_back(std::stod(item));
                        }
                        break;
                }
            }

            // find where the start of the scheme block is
            auto start_pos = std::ranges::search(line, std::format("BEGIN {}", scheme));

            if (!start_pos.empty()) {
                start_line = i;
            }
            i+=1;
        }

        if (debug){
            std::println("atomic_numbers: {}", iso.atomic_numbers);
            std::println("elements: {}",iso.elements);
            std::println("mass_numbers: {}", iso.mass_numbers);
            std::println("percentages: {}", iso.percentages);
        }

        return iso;
    }
}

namespace fourdst::composition {
    Composition get_composition_record(const std::string& metal_fraction_scheme,
                                                                        const std::string& isotopic_percentage_scheme,
                                                                        double initial_z, double initial_y) {


        std::vector<char> data;

        io::ChemicalFileParser parser;
        io::CompositionData metals;

        io::IsotopicPercentage isotopes;

        data = std::ranges::to<std::vector<char>>(StandardMetalFractions);

        metals = parser.parse_compositon_data(data,metal_fraction_scheme);
        isotopes = parser.parse_isotopic_percentage(data,isotopic_percentage_scheme);

        std::string name;
        std::vector<atomic::Species> species;


        // construct name of the isotopes for all elements
        for (const auto [E,A] : std::ranges::views::zip(isotopes.elements, isotopes.mass_numbers)){
            if (std::ranges::contains(metals.elements,E ) || E == "H" || E == "He") {
                name = std::format("{}-{}",E,A);
                auto SpeciesObject = atomic::species.at(name);
                species.push_back(SpeciesObject);
                // std::println("Species: {} has mass: {}", SpeciesObject.name(), SpeciesObject.mass());
            }
        }

        std::vector<double> massFracs;
        std::unordered_map<std::string, double> metal_fractions;

        // hydrogen and helium are treated separately
        // H1
        massFracs.push_back(std::max(0.0, std::min(1.0, 1.0 - (initial_z + initial_y))));
        // H2
        massFracs.push_back(0.0);
        // He3
        // anders & grevesse 1989 solar mass fractions
        double xsol_he3,xsol_he4;
        xsol_he3 = 2.9291e-05;
        xsol_he4 = 2.7521e-01;
        massFracs.push_back(initial_y*xsol_he3/(xsol_he3 + xsol_he4));
        // He4
        massFracs.push_back(initial_y*xsol_he4/(xsol_he3 + xsol_he4));
        // Metals

        double ztotal = 1.0-std::accumulate(massFracs.begin(), massFracs.end(), 0.0);

        // multiply by atomic weight if needed

        if (metals.requires_atomic_weight){
            // get isotope with max abundance for each metal
            // and store the corresponding mass number
            auto element_atomic_weight = [&isotopes]() {
                std::unordered_map<std::string, std::pair<double, int>> elem_info;

                for (const auto& [iso, prcnt, a] : std::views::zip(isotopes.elements, isotopes.percentages, isotopes.mass_numbers)) {
                    if (iso == "H" || iso == "He") {
                        continue;
                    }
                    if (elem_info.contains(iso) && elem_info.at(iso).first <= prcnt) {
                        elem_info[iso] = std::make_pair(prcnt, a);
                    } else if (! elem_info.contains(iso)) {
                        elem_info[iso] = std::make_pair(prcnt, a);
                    }
                }
                return elem_info;
            }();

            for (const auto [E,A] : std::ranges::views::zip(metals.elements, metals.abundances)) {
                // std::println("element: {}", E);
                auto name = std::format("{}-{}",E,element_atomic_weight.at(E).second);
                // std::println("{}", name);
                auto SpeciesObject = atomic::species.at(name);
                double weight = SpeciesObject.mass();
                metal_fractions.emplace(E,A*weight);
                // std::println("End");
            }
        } else {
            for (const auto [E,A] : std::ranges::views::zip(metals.elements, metals.abundances)) {
                metal_fractions.emplace(E,A);
            }
        }

        double sum = [&metal_fractions]() {
          double accumulator = 0.0;
            for (const auto& frac : metal_fractions | std::views::values) {
                accumulator += frac;
            }
            return accumulator;
        }();

        for (auto& frac : metal_fractions | std::views::values) {
            frac/=sum;
        }

        double zsum = 0.0;

        // get mass Fracs for each metal and scale it to required ztotal
        for (size_t i = 0; i < species.size();++i) {
            size_t Z = isotopes.atomic_numbers[i];
            if (Z<=2) continue;

            if (metal_fractions.contains(isotopes.elements[i])) {
                double frac = 1e-2*isotopes.percentages[i]*species[i].mass();
                double frac_sum = 0.0;
                size_t j = i ;
                while (j<species.size() && isotopes.atomic_numbers[j] == Z) {
                    frac_sum += 1e-2*isotopes.percentages[j]*species[j].mass();
                    ++j;
                }
                // extract zfrac for the corresponding Z symbol/ isotopes.elements
                double zfrac = metal_fractions.at(isotopes.elements[i]);
                auto temp = ztotal*zfrac*frac/frac_sum;
                massFracs.push_back(temp);
                zsum += temp;
                // std::println("isotope:{}",species[i].name());
            }
        }

        // std::println("ztotal: {}, zsum:{}", ztotal, zsum);
        //Renormalize
        if (zsum > 0.0) {
            for (size_t i = 0; i < massFracs.size();++i) {
                if (isotopes.atomic_numbers[i]<=2) continue;
                massFracs[i] *= ztotal/zsum;
            }
        }

        //
        Composition comp = buildCompositionFromMassFractions(species, massFracs);
        return comp;

    }

    Composition get_composition_record(const io::SolarCompositions metal_fraction_scheme,
                                       const io::IsotopicPercentages isotopic_percentage_scheme,
                                       double initial_z,
                                       double initial_y) {
        return get_composition_record(
            io::SolarCompositions_to_string_map.at(metal_fraction_scheme),
            io::IsotopicPercentages_to_string_map.at(isotopic_percentage_scheme),
            initial_z,
            initial_y
        );
    }
}
