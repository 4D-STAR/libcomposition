![logo](assets/logo/logo.png)
[documentation](https://4d-star.github.io/libcomposition/html/)

# Introduction

`libcomposition` is a modern, C++23 library, for the creation, manipulation, and analysis of astrophysical chemical
compositions. It provides a robust and type‑safe interface for assembling a set of isotopes together with their molar
abundances and for deriving commonly used bulk properties (mass fractions, number fractions, canonical X/Y/Z, mean
particle mass, and electron abundance). `libcomposition` is designed to be tighly integrated into SERiF and related
projects such as GridFire.

### Key Features

- **Type–Safe Species Representation**: Strongly typed isotopes (`fourdst::atomic::Species`) generated from evaluated nuclear data (AME2020 / NUBASE2020).
- **Molar Abundance Core**: Stores absolute molar abundances and derives all secondary quantities (mass / number fractions, mean particle mass, electron abundance) on demand, with internal caching.
- **Canonical Composition Support**: Direct computation of canonical (X: Hydrogen, Y: Helium, Z: Metals) mass fractions via `getCanonicalComposition()`.
- **Convenience Construction**: Helper utilities for constructing compositions from a vector or set of mass fractions (`buildCompositionFromMassFractions`).
- **Deterministic Ordering**: Species are always stored and iterated lightest→heaviest (ordering defined by atomic mass) enabling uniform vector interfaces.
- **Clear Exception Hierarchy**: Explicit error signaling for invalid symbols, unregistered species, and inconsistent input data.
- **Meson + pkg-config Integration**: Simple build, install, and consumption in external projects.

---

# Installation

libcomposition can be installed either from source or as part of the `fourdst` project.

`libcomposition` uses the Meson build system. A C++23 compatible compiler is required.

### Build Steps

**Setup the build directory:**

The first step is to use meson to set up an out of source build. Note that this means that you can have multiple builds configured and cleanly separated!

```bash
meson setup builddir
```

**Compile the library:**

meson by default uses ninja to compile so it should be very fast; however, gcc is very slow when compiling the species database so that might take some time (clang tends to be very fast for this).

```bash
meson compile -C builddir
```

**Install the library:**

This will also install a pkg-config file!

```bash
sudo meson install -C builddir
```

### Build Options

You can enable the generation of a `pkg-config` file during the setup step, which simplifies linking the library in other projects. By default this is true; it can be useful to disable this when using some build system orchestrator (such as meson-python).

```bash
# Enable pkg-config file generation
meson setup builddir -Dpkg-config=true
```

---

# Usage

Below are focused examples illustrating the current API. All examples assume headers are available via pkg-config or your include path.

#### 1. Constructing a Composition from Symbols

```cpp
#include <iostream>
#include "fourdst/composition/composition.h"

int main() {
    using namespace fourdst::composition;

    // Register symbols upon construction (no molar abundances yet -> default 0.0)
    Composition comp({"H-1", "He-4", "C-12"});

    // Set molar abundances (absolute counts; they need not sum to 1.0)
    comp.setMolarAbundance("H-1", 10.0);
    comp.setMolarAbundance("He-4", 3.0);
    comp.setMolarAbundance("C-12", 0.25);

    // Query derived properties
    double x_h1 = comp.getMassFraction("H-1");
    double y_he4 = comp.getNumberFraction("He-4");
    auto canon = comp.getCanonicalComposition(); // X, Y, Z mass fractions

    std::cout << "H-1 mass fraction: " << x_h1 << "\n";
    std::cout << "He-4 number fraction: " << y_he4 << "\n";
    std::cout << canon << "\n"; // <CanonicalComposition: X=..., Y=..., Z=...>
}
```

#### 2. Constructing from Strongly Typed Species

```cpp
#include <iostream>
#include "fourdst/composition/composition.h"
#include "fourdst/atomic/species.h"

int main() {
    using namespace fourdst::composition;
    using namespace fourdst::atomic;

    // Build directly from species constants
    Composition comp(std::vector<Species>{H_1, He_4, O_16});

    comp.setMolarAbundance(H_1, 5.0);
    comp.setMolarAbundance(He_4, 2.5);
    comp.setMolarAbundance(O_16, 0.1);

    std::cout << "Mean particle mass: " << comp.getMeanParticleMass() << " g/mol\n";
    std::cout << "Electron abundance (Ye): " << comp.getElectronAbundance() << "\n";
}
```

#### 3. Building from Mass Fractions (Helper Utility)

```cpp
#include <iostream>
#include "fourdst/composition/utils.h"

int main() {
    using namespace fourdst::composition;

    std::vector<std::string> symbols = {"H-1", "He-4", "C-12"};
    std::vector<double> mf       = {0.70, 0.28, 0.02}; // Must sum to ~1 within tolerance

    Composition comp = buildCompositionFromMassFractions(symbols, mf);

    auto canon = comp.getCanonicalComposition();
    std::cout << canon << "\n";
}
```

#### 4. Iterating and Sorted Vector Interfaces

```cpp
#include <iostream>
#include "fourdst/composition/composition.h"

int main() {
    using namespace fourdst::composition;

    Composition comp({"H-1", "C-12", "He-4"}); // Internally sorted by mass (H < He < C)
    comp.setMolarAbundance({"H-1", "He-4", "C-12"}, {10.0, 3.0, 0.25});

    // Ordered iteration (lightest -> heaviest)
    for (const auto &[sp, y] : comp) {
        std::cout << sp << ": molar = " << y << "\n";
    }

    // Vector access (index corresponds to ordering by atomic mass)
    auto molarVec = comp.getMolarAbundanceVector();
    auto massVec  = comp.getMassFractionVector();

    size_t idx_he4 = comp.getSpeciesIndex("He-4");
    std::cout << "He-4 index: " << idx_he4 << ", molar abundance at index: " << molarVec[idx_he4] << "\n";
}
```

#### 5. Accessing Specific Derived Quantities

```cpp
// Assume 'comp' is already populated.

double mf_c12   = comp.getMassFraction("C-12");
double nf_c12   = comp.getNumberFraction("C-12");
double mol_c12  = comp.getMolarAbundance("C-12");
double meanA    = comp.getMeanParticleMass();
double Ye       = comp.getElectronAbundance();
auto   canon    = comp.getCanonicalComposition();
```

#### 6. Exception Handling Examples

```cpp
#include <iostream>
#include "fourdst/composition/composition.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"

int main() {
    using namespace fourdst::composition;
    using namespace fourdst::composition::exceptions;

    Composition comp;

    try {
        // Unknown symbol (not in species database)
        comp.registerSymbol("Xx-999");
    } catch (const UnknownSymbolError &e) {
        std::cerr << "Caught UnknownSymbolError: " << e.what() << "\n";
    }

    comp.registerSymbol("H-1");
    try {
        // Unregistered symbol used in a setter
        comp.setMolarAbundance("He-4", 1.0); // He-4 not registered yet
    } catch (const UnregisteredSymbolError &e) {
        std::cerr << "Caught UnregisteredSymbolError: " << e.what() << "\n";
    }

    comp.registerSymbol("He-4");
    try {
        comp.setMolarAbundance("H-1", -3.0);
    } catch (const InvalidCompositionError &e) { 
        std::cerr << "Caught InvalidCompositionError: " << e.what() << "\n";
    }

    // Mass fraction construction validation
    try {
        Composition bad = buildCompositionFromMassFractions({"H-1", "He-4"}, {0.6, 0.5}); // sums to 1.1
    } catch (const InvalidCompositionError &e) {
        std::cerr << "Caught InvalidCompositionError: " << e.what() << "\n";
    }
}
```

---

@section exceptions_sec Possible Exception States

The library surfaces errors through a focused hierarchy in `fourdst::composition::exceptions`:

| Exception Type | When It Occurs |
|----------------|----------------|
| `UnknownSymbolError` | A string symbol does not correspond to any known isotope in the compiled species database. |
| `UnregisteredSymbolError` | A valid species/symbol is used before being registered with a Composition instance. |
| `InvalidCompositionError` | Construction from mass fractions fails validation (sum deviates from unity beyond tolerance) or canonical (X+Y+Z) check fails. |
| `CompositionError` | Base class; may be thrown for generic composition-level issues (e.g. negative abundances via the documented `InvalidAbundanceError` contract). |

Recommended patterns:
- Validate externally provided symbol lists before calling bulk registration.
- Use species‑based overloads (strongly typed) where possible for slightly lower overhead (no symbol resolution).
- Wrap construction from mass fractions in a try/catch to surface normalization issues early.

---

# Linking and Integration

### Linking with pkg-config

If you installed `libcomposition` with the `pkg-config` option enabled, you can get the necessary compiler and linker flags easily:

```bash
# Get compiler flags (include paths)
pkg-config --cflags fourdst_composition

# Get linker flags (library paths and names)
pkg-config --libs fourdst_composition
```

**Example compilation command:**
```bash
g++ my_app.cpp $(pkg-config --cflags --libs fourdst_composition) -o my_app
```