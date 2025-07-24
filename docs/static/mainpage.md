@mainpage libcomposition: A Modern C++ Library for Chemical Compositions

@section intro_sec Introduction

`libcomposition` is a modern C++23 library designed for the creation, manipulation, and analysis of chemical compositions, with a focus on astrophysical applications. It provides a robust and user-friendly interface for handling material compositions defined by mass or number fractions.

### Key Features

-   **Dual-Mode Operation**: Natively supports compositions defined by **mass fraction** or **number fraction**.
-   **Rich Atomic Database**: Includes a comprehensive, header-only database of isotopic properties (mass, half-life, spin, etc.) generated from the AME2020 and NUBASE2020 evaluations.
-   **Type Safety and Error Handling**: Utilizes a clear exception hierarchy to report errors, such as using an unregistered isotope or accessing data from a non-validated composition.
-   **Powerful Functionality**: Core features include mixing, subsetting, and on-the-fly conversion between mass and number fractions.
-   **Easy Integration**: Designed for seamless integration with other projects using the Meson build system and `pkg-config`.

---

@section install_sec Installation

`libcomposition` uses the Meson build system. A C++23 compatible compiler is required.

### Build Steps

**Setup the build directory:**

The first step is to use meson to set up an out of source build. Note that this means that you can have multiple builds configured and cleanly seperated!

```bash
meson setup builddir
```

**Compile the library:**

meson by default uses ninja to compile so it should be very fast; however, gcc is very slow when compiling the species database so that migth take some time (clang tends to be very fast for this).

```bash
meson compile -C builddir
```

 **Install the library:**

This will also install a pkg-config file!

```bash
sudo meson install -C builddir
```

### Build Options

You can enable the generation of a `pkg-config` file during the setup step, which simplifies linking the library in other projects. by default this is true; it can be useful to disable this when using some build system orgestrator (such as meson-python).

```bash
# Enable pkg-config file generation
meson setup builddir -Dpkg-config=true
```

---

@section usage_sec Usage

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

### C++ Usage Examples

#### 1. Basic Mass Fraction Composition

The most common use case is defining a composition by mass fractions (X, Y, Z).

```cpp
#include <iostream>
#include "fourdst/composition/composition.h"

int main() {
    // 1. Create a composition object
    fourdst::composition::Composition comp;

    // 2. Register the symbols you want to use
    comp.registerSymbol("H-1");
    comp.registerSymbol("He-4");

    // 3. Set their mass fractions
    comp.setMassFraction("H-1", 0.75);
    comp.setMassFraction("He-4", 0.25);

    // 4. Finalize the composition to validate it and compute global properties
    if (comp.finalize()) {
        std::cout << "Composition finalized successfully!" << std::endl;
        std::cout << "H-1 Mass Fraction: " << comp.getMassFraction("H-1") << std::endl;
        std::cout << "Mean Particle Mass: " << comp.getMeanParticleMass() << " g/mol" << std::endl;
    } else {
        std::cerr << "Failed to finalize composition." << std::endl;
    }

    return 0;
}
```

#### 2. Number Fraction Composition and Mode Switching

The library can also work with number (mole) fractions and switch between modes.

```cpp
#include "fourdst/composition/composition.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"

void number_fraction_example() {
    fourdst::composition::Composition comp;

    // Register symbols in number fraction mode
    comp.registerSymbol("H-1", false); // massFracMode = false
    comp.registerSymbol("He-4", false);

    comp.setNumberFraction("H-1", 0.9);
    comp.setNumberFraction("He-4", 0.1);

    if (comp.finalize()) {
        // We can get number fractions directly
        std::cout << "He-4 Number Fraction: " << comp.getNumberFraction("He-4") << std::endl;

        // Or get the equivalent mass fraction
        std::cout << "He-4 Mass Fraction: " << comp.getMassFraction("He-4") << std::endl;

        // Switch the entire composition to mass fraction mode
        comp.setCompositionMode(true); // true for mass fraction mode

        // Now, getting the mass fraction is a direct lookup
        std::cout << "He-4 Mass Fraction (after mode switch): " << comp.getMassFraction("He-4") << std::endl;
    }
}
```

#### 3. Mixing Two Compositions

You can easily mix two compositions. The library handles the union of all species.

```cpp
#include "fourdst/composition/composition.h"

void mixing_example() {
    // Composition 1: Pure Hydrogen
    fourdst::composition::Composition comp1({"H-1"}, {1.0});

    // Composition 2: Pure Helium
    fourdst::composition::Composition comp2({"He-4"}, {1.0});

    // Mix them with a 50/50 ratio using the '+' operator
    fourdst::composition::Composition mixed = comp1 + comp2;

    // Mix them with a 75/25 ratio using the mix() method
    // 0.75 of comp1, 0.25 of comp2
    fourdst::composition::Composition mixed2 = comp1.mix(comp2, 0.75);

    std::cout << "50/50 Mix H-1: " << mixed.getMassFraction("H-1") << std::endl;   // -> 0.5
    std::cout << "75/25 Mix H-1: " << mixed2.getMassFraction("H-1") << std::endl;  // -> 0.75
}
```

#### 4. Error Handling

The library uses exceptions to report errors. Always wrap calls in a `try-catch` block for robust code.

```cpp
#include "fourdst/composition/composition.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"

void error_example() {
    fourdst::composition::Composition comp;
    comp.registerSymbol("H-1");
    comp.setMassFraction("H-1", 1.0);

    try {
        // This will throw, because the composition is not finalized yet.
        double mass = comp.getMassFraction("H-1");
    } catch (const fourdst::composition::exceptions::CompositionNotFinalizedError& e) {
        std::cerr << "Caught expected error: " << e.what() << std::endl;
    }

    try {
        // This will throw, because "Li-6" was never registered.
        comp.setMassFraction("Li-6", 0.1);
    } catch (const fourdst::composition::exceptions::UnregisteredSymbolError& e) {
        std::cerr << "Caught expected error: " << e.what() << std::endl;
    }
}
```

#### 5. Accessing Atomic Data

You can directly access the static database of all known species.

```cpp
#include "fourdst/composition/species.h" // Provides static instances like H_1
#include "fourdst/composition/atomicSpecies.h" // Provides the main 'species' map

void data_example() {
    // Access via the map
    const auto& fe56 = fourdst::atomic::species.at("Fe-56");
    std::cout << "Fe-56 mass: " << fe56->mass() << std::endl;

    // Access via the static instance
    std::cout << "H-1 spin: " << fourdst::atomic::H_1.spin() << std::endl;
    std::cout << "F-18 half-life: " << fourdst::atomic::F_18.halfLife() << " s" << std::endl;
}
```

---

@section test_sec Testing

`libcomposition` is tested using the GoogleTest framework. The test suite provides high coverage of the library's functionality.

### Test Coverage Includes:

-   **Atomic Data Validation**: Spot checks on isotopic properties (mass, half-life, spin) for a wide range of elements to ensure the underlying data files are parsed and represented correctly.
-   **Core `Composition` Workflow**: Verification of object construction, symbol registration (for both valid and invalid symbols), and the complete workflow of setting and getting both mass and number fractions.
-   **Finalization Logic**: Ensures that `finalize()` is a required step before querying data. Tests the validation logic for compositions that sum to 1.0 and the auto-normalization feature (`finalize(true)`).
-   **Advanced Features**: Dedicated tests for `mix()`, `subset()`, `setCompositionMode()`, and the calculation of derived quantities like `getMolarAbundance()` and `getMeanAtomicNumber()`.
-   **Exception Handling**: Confirms that invalid operations (e.g., using an unregistered symbol, mixing un-finalized compositions) correctly throw exceptions from the `fourdst::composition::exceptions` hierarchy.

---

@section api_sec API Reference

For a complete list of all classes, methods, and functions, please see the **<a href="namespaces.html">Namespaces</a>** and **<a href="annotated.html">Classes</a>** sections of this documentation.
