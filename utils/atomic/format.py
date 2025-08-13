import pandas as pd
import numpy as np
import argparse
import os


# Column specifications for the NUBASE2020 data file.
# These are character positions for fixed-width fields.
nubase_col_specs = [
    (0, 3),    # Mass number (A)
    (4, 8),    # ZZZi identifier for isomer level
    (11, 16),  # A_El, e.g., "56Fe"
    (69, 78),  # Half-life value
    (79, 81),  # Half-life unit
    (89, 102), # Spin and parity
    (120, 209) # Decay modes
]
nubase_column_names = [
    "a",             
    "ZZZi",
    "A_El",          
    "halfLife",
    "halfLifeUnit",
    "spinParity",
    "decayModes"
]

# Column specifications for the AME2020 data file.
ame_col_specs = [
    (0, 1),
    (1, 4),
    (4, 9),
    (9, 14),
    (14, 19),
    (20, 23),
    (23, 27),
    (28, 42),
    (42, 54),
    (54, 67),
    (68, 78),
    (79, 81),
    (81, 94),
    (94, 105),
    (106, 109),
    (110, 123),
    (123, 135)
]
ame_column_names = [
    "control",
    "nz",
    "n",
    "z",
    "a",
    "el",
    "o",
    "massExcess",
    "massExcessUnc",
    "bindingEnergy",
    "bindingEnergyUnc",
    "betaCode",
    "betaDecayEnergy",
    "betaDecayEnergyUnc",
    "atomicMassInt",
    "atomicMassFrac",
    "atomicMassUnc"
]


def combine_atomic_mass(row):
    """
    Combines the integer and fractional parts of the atomic mass from an AME2020 row.

    The AME2020 format provides atomic mass as two separate columns: an integer part
    and a fractional part in micro-u. This function merges them into a single
    double-precision float.

    Args:
        row (pd.Series): A row from the AME DataFrame.

    Returns:
        float or None: The combined atomic mass, or None if conversion fails.
    """
    intPart = str(row["atomicMassInt"]).strip()
    fracPart = str(row["atomicMassFrac"]).strip()
    try:
        # The fractional part is in micro-u, so it's divided by 1e6.
        combined = int(intPart) + float(fracPart)/1e6
        return combined
    except ValueError:
        return None

def mkInstanceName(row):
    """
    Creates a valid C++ variable name for a species.

    Args:
        row (pd.Series): A row from the DataFrame.

    Returns:
        str: A C++-safe instance name (e.g., "Fe_56").
    """
    speciesName = f"{row['el'].strip()}-{row['a']}"
    return speciesName.replace("-", "_")

def convert_half_life_to_seconds(row):
    """
    Converts the half-life value and unit from a NUBASE2020 row into seconds.

    Handles various time units (from yoctoseconds to petayears) and the 'stbl'
    (stable) identifier.

    Args:
        row (pd.Series): A row from the NUBASE DataFrame.

    Returns:
        float: The half-life in seconds. Returns np.inf for stable isotopes
               and 0.0 for parsing errors.
    """
    value_str = str(row['halfLife']).strip()
    unit = str(row['halfLifeUnit']).strip()

    if value_str == 'stbl':
        return np.inf

    try:
        # The value string can contain non-numeric characters like '#' for uncertainty.
        value = float(value_str.replace('#', ''))
    except ValueError:
        return 0.0 

    # Conversion factors from various units to seconds.
    factors = {
        'ys': 1e-24, 'zs': 1e-21, 'as': 1e-18, 'fs': 1e-15,
        'ps': 1e-12, 'ns': 1e-9,  'us': 1e-6,  'ms': 1e-3,
        's': 1.0, 'm': 60.0, 'h': 3600.0, 'd': 86400.0,
        'y': 3.15576e7,  
        'ky': 3.15576e10, 'My': 3.15576e13, 'Gy': 3.15576e16,
        'Ty': 3.15576e19, 'Py': 3.15576e22, 'Ey': 3.15576e25
    }
    return value * factors.get(unit, 0.0)


def formatSpecies(row):
    """
    Formats a DataFrame row into a C++ static Species object instantiation string.

    Args:
        row (pd.Series): A row from the merged DataFrame.

    Returns:
        tuple[str, str]: A tuple containing the full C++ instantiation string
                         and the instance name.
    """
    name = f"{row['el'].strip()}-{row['a']}"
    instanceName = name.replace("-", "_")
    nz = int(row['nz'])
    n = int(row['n'])
    z = int(row['z'])
    a = int(row['a'])
    bindingEnergy = float(row['bindingEnergy'])
    atomicMass = float(row['atomicMass'])
    atomicMassUnc = float(row['atomicMassUnc'])
    NaN = "std::numeric_limits<double>::quiet_NaN()"
    
    try:
        # Clean up beta decay energy string before converting to float.
        betaDecayEnergy = float(row['betaDecayEnergy'].replace("#", "").replace("*", ""))
    except (ValueError, AttributeError):
        betaDecayEnergy = NaN

    ### --- Add new NUBASE fields ---
    halfLife_s = row.get('halfLife_s', 'std::numeric_limits<double>::infinity()')
    if halfLife_s == np.inf:
        halfLife_s = 'std::numeric_limits<double>::infinity()'
        
    # Escape double quotes in string fields for C++ compatibility.
    spinParity = str(row.get('spinParity', '')).strip().replace('"', '\\"')
    decayModes = str(row.get('decayModes', '')).strip().replace('"', '\\"')

    instantiation = (
        f"static const Species {instanceName}(\"{name}\", \"{row['el'].strip()}\", {nz}, {n}, {z}, {a}, "
        f"{bindingEnergy}, \"{row['betaCode']}\", {betaDecayEnergy}, "
        f"{halfLife_s}, \"{spinParity}\", \"{decayModes}\", "
        f"{atomicMass}, {atomicMassUnc});"
    )
    return instantiation, instanceName

def formatSpeciesDefines(row):
    """
    Generates C++ preprocessor define guards for a species.
    Note: This function is not currently used in the script's main execution path.

    Args:
        row (pd.Series): A row from the DataFrame.

    Returns:
        str: A string containing the #ifndef/#define block.
    """
    instanceName = f"SERIF_SPECIES_{formatSpecies(row)[1]}"
    define = f"""#ifndef {instanceName.upper()}
    #define {instanceName.upper()}
#endif // {instanceName.upper()}"""
    return define


def formatHeader(dataFrame):
    """
    Generates the complete C++ header file content as a string.

    This function creates all the static Species instantiations and builds an
    unordered_map to provide string-based lookup at runtime.

    Args:
        dataFrame (pd.DataFrame): The final merged DataFrame containing all species data.

    Returns:
        str: The content of the C++ header file.
    """
    header = f"""#pragma once
#include <unordered_map>
#include <string_view>
#include <string>
#include <limits> // Required for std::numeric_limits
#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/elements.h"

namespace fourdst::atomic {{
    // Instantiate all species as static const objects.
    {'\n    '.join([formatSpecies(row)[0] for index, row in dataFrame.iterrows()])}
    
    // Create a map from species name (e.g., "H-1") to a pointer to the species object.
    static const std::unordered_map<std::string, const Species*> species = {{
        {'\n        '.join([f'{{"{row["el"].strip()}-{row["a"]}", {mkInstanceName(row)}}},' for index, row in dataFrame.iterrows()])}
    }};
    Species az_to_species(const int a, const int z) {{
    const std::string element_symbol = element_symbol_map.at(static_cast<uint8_t>(z));
        const std::string species_symbol = element_symbol + "-" + std::to_string(a);
        return species.at(species_symbol);
    }};

}}; // namespace fourdst::atomic
"""
    return header

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert AME2020 and NUBASE2020 data to a C++ header file.")
    parser.add_argument("ame_input", help="Input file path for AME2020 (mass.mas20).")
    parser.add_argument("nubase_input", help="Input file path for NUBASE2020.")
    parser.add_argument("-o", "--output", help="Output file path.", default="species.h")
    args = parser.parse_args()

    for path in [args.ame_input, args.nubase_input]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

    # Read the AME2020 mass data, skipping the header.
    ame_df = pd.read_fwf(args.ame_input, colspecs=ame_col_specs, names=ame_column_names, skiprows=36, keep_default_na=False)
    ame_df["atomicMass"] = ame_df.apply(combine_atomic_mass, axis=1)
    ame_df.drop(columns=["atomicMassInt", "atomicMassFrac"], inplace=True)
    ame_df['el'] = ame_df['el'].str.strip()

    # Read the NUBASE2020 nuclear properties data, skipping the header.
    nubase_df = pd.read_fwf(args.nubase_input, colspecs=nubase_col_specs, names=nubase_column_names, skiprows=1, keep_default_na=False)
    nubase_df['a'] = pd.to_numeric(nubase_df['a'], errors='coerce')
    nubase_df.dropna(subset=['a'], inplace=True)
    nubase_df['a'] = nubase_df['a'].astype(int)
    # Extract the element symbol (e.g., 'Fe' from '56Fe').
    nubase_df['el'] = nubase_df['A_El'].str.extract(r'([a-zA-Z]+)')[0].str.strip()
    # Filter for ground states only. The 4th char in ZZZi is ' ' or '0' for ground states.
    nubase_df = nubase_df[nubase_df['ZZZi'].str[3:4].isin([' ', '0'])]
    # In case of multiple entries for the same isotope (e.g., different isomers), keep the first one.
    nubase_df.drop_duplicates(subset=['a', 'el'], keep='first', inplace=True)
    nubase_df['halfLife_s'] = nubase_df.apply(convert_half_life_to_seconds, axis=1)
    

    print("--- AME DataFrame ---")
    print(ame_df[['a', 'el']].head())
    print("\n--- NUBASE DataFrame ---")
    print(nubase_df[['a', 'el', 'halfLife_s']].head())
    print("\n")

    # Merge the AME (mass) and NUBASE (properties) data.
    # A left merge keeps all isotopes from AME and adds properties from NUBASE where available.
    merged_df = pd.merge(ame_df, nubase_df[['a', 'el', 'halfLife_s', 'spinParity', 'decayModes']], on=['a', 'el'], how='left')

    print("--- Merged DataFrame ---")
    print(merged_df[['a', 'el', 'halfLife_s']].head(10))
    print("\n--- Merge Stats ---")
    print(f"Total rows in final data: {len(merged_df)}")
    print(f"Number of rows with a valid half-life after merge: {merged_df['halfLife_s'].notna().sum()}")
    print("-" * 20)

    # Generate the C++ header content from the merged data.
    header = formatHeader(merged_df)
    with open(args.output, "w") as f:
        f.write(header)

    print(f"Successfully generated C++ header at {args.output}")
