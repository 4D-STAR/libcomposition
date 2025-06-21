import pandas as pd

# Define fixed-width column specifications based on the format:
# a1 (width 1), i3 (width 3), i5 (width 5), i5 (width 5), i5 (width 5),
# 1x (skip 1), a3 (width 3), a4 (width 4), 1x (skip 1),
# f14.6 (width 14), f12.6 (width 12), f13.5 (width 13),
# 1x (skip 1), f10.5 (width 10), 1x (skip 1),
# a2 (width 2), f13.5 (width 13), f11.5 (width 11),
# 1x (skip 1), i3 (width 3), 1x (skip 1),
# f13.6 (width 13), f12.6 (width 12)
# Compute cumulative positions (0-indexed):
colSpecs = [
    (0, 1),    # control
    (1, 4),    # NZ
    (4, 9),    # N
    (9, 14),   # Z
    (14, 19),  # A
    # skip 1 char at position 19; next field starts at 20
    (20, 23),  # el
    (23, 27),  # o
    # skip 1 char at position 27; next field starts at 28
    (28, 42),  # massExcess (f14.6)
    (42, 54),  # massExcessUnc (f12.6)
    (54, 67),  # bindingEnergy (f13.5)
    # skip 1 char at position 67; next field starts at 68
    (68, 78),  # bindingEnergyUnc (f10.5)
    # skip 1 char at position 78; next field starts at 79
    (79, 81),  # betaCode (a2)
    (81, 94),  # betaDecayEnergy (f13.5)
    (94, 105), # betaDecayEnergyUnc (f11.5)
    # skip 1 char at position 105; next field starts at 106
    (106, 109),# atomicMassInt (i3)
    # skip 1 char at position 109; next field starts at 110
    (110, 123),# atomicMassFrac (f13.6)
    (123, 135) # atomicMassUnc (f12.6)
]

# Define column names (using camelCase for variables)
columnNames = [
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
    Combine the integer and fractional parts of the atomic mass.
    For example, if atomicMassInt is '1' and atomicMassFrac is '008664.91590',
    this function returns float('1008664.91590').
    """
    intPart = str(row["atomicMassInt"]).strip()
    fracPart = str(row["atomicMassFrac"]).strip()
    try:
        combined = int(intPart) + float(fracPart)/1e6
        return combined
    except ValueError:
        return None

def mkInstanceName(row):
    """
    Make a c++ instance name from the element and atomic number.
    """
    speciesName = f"{row['el'].strip()}-{row['a']}"
    return speciesName.replace("-", "_")

def formatSpecies(row):
    """
    Format c++ instantiation of Species struct from row data.
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
        betaDecayEnergy = float(row['betaDecayEnergy'].replace("#", "").replace("*", ""))
    except ValueError:
        betaDecayEnergy = NaN
    instantiation = f"static const Species {instanceName}(\"{name}\", \"{row['el']}\", {nz}, {n}, {z}, {a}, {bindingEnergy}, \"{row['betaCode']}\", {betaDecayEnergy}, {atomicMass}, {atomicMassUnc});"
    return instantiation, instanceName

def formatSpeciesDefines(row):
    instanceName = f"SERIF_SPECIES_{formatSpecies(row)[1]}"
    define = f"""#ifndef {instanceName.upper()}
    #define {instanceName.upper()}
#endif // {instanceName.upper()}"""
    return define

def formatHeader(dataFrame):
    """
    Format c++ header file from DataFrame.
    """
    header = f"""#pragma once
#include <unordered_map>
#include <string_view>
#include <string>
#include "atomicSpecies.h"

namespace fourdst::atomic {{
    {'\n    '.join([formatSpecies(row)[0] for index, row in dataFrame.iterrows()])}
    static const std::unordered_map<std::string, Species> species = {{
        {'\n        '.join([f'{{"{row["el"].strip()}-{row["a"]}", {mkInstanceName(row)}}},' for index, row in dataFrame.iterrows()])}
    }};
}}; // namespace fourdst::atomic

{'\n'.join([formatSpeciesDefines(row) for index, row in dataFrame.iterrows()])}
"""
    return header

if __name__ == "__main__":
    import argparse 
    import os
    parser = argparse.ArgumentParser(description="Convert mass data to c++ header file.")
    parser.add_argument("input", help="Input file path.")
    parser.add_argument("-o", "--output", help="Output file path.", default="../../assets/static/atomic/include/species.h")
    args = parser.parse_args()


    if not os.path.exists(args.input):
        raise FileNotFoundError(f"File not found: {args.input}")

    # Read the file (adjust the skiprows value if your header differs)
    dataFrame = pd.read_fwf(args.input, colspecs=colSpecs, names=columnNames, skiprows=36)

    # Combine the two atomic mass fields into one float column
    dataFrame["atomicMass"] = dataFrame.apply(combine_atomic_mass, axis=1)
    dataFrame.drop(columns=["atomicMassInt", "atomicMassFrac"], inplace=True)

    # Format the header
    header = formatHeader(dataFrame)
    with open(args.output, "w") as f:
        f.write(header)
