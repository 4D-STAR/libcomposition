## Unreleased

## v2.0.0 (2025-11-07)

### BREAKING CHANGE

- The entire old API has been broken. There is no longer any need to finalize. In fact the entire concept of finalization has been removed. Further, the entire CompositionEntry and GlobalComposition data structure has been removed. Any code written for the old version will no longer work and major reworking will be needed to use the new version.

### Feat

- **composition**: changed how composition is conmstructed

## v1.9.1 (2025-10-22)

### Perf

- **Composition**: improved static correctness
- **Composition**: changed logger aquisition to static to prevent each Composition object from needing to instantiate a logger whenever its built
- **CompositionEntry**: made symbol and isotope std::optional for comopsition entry so that construction can be cheaper

## v1.9.0 (2025-10-12)

### Feat

- **Composition**: Composition now inherits from abstract base class

## v1.8.2 (2025-10-09)

### Perf

- **az_to_species**: marked noexcept

## v1.8.1 (2025-10-08)

### Feat

- **az_to_species**: az_to_species now returns an expected and error type

## v1.8.0 (2025-10-06)

### Feat

- **meson.build**: version bump (1.7.0 -> 1.8.0)
- **Added-ability-to-get-electron-abundance-and-fixed-some-conversion-bugs**: Now Ye can be retrived directly from the composition object. Further a bug which prevented proper conversion to and from number or mass frac modes without messing up the numbers has been resolved

## v1.7.0 (2025-09-16)

### Feat

- **composition**: added uniform tools to get vector representation of mass fraction, number fraction, and molar abundance

## v1.6.0 (2025-08-13)

### Feat

- **species-lookup**: added function to get species from a and z

## v1.5.2 (2025-07-24)

## v1.5.1 (2025-07-22)

## v1.5.0 (2025-07-22)

## v1.4.1 (2025-07-21)

### Feat

- **composition**: added more expressive errors

## v1.4.0 (2025-07-14)

### Feat

- **composition**: added species queries and < > operators for species based on mass

## v1.3.0 (2025-07-04)

### Feat

- **composition**: added stl compatible iterator

## v1.2.0 (2025-07-02)

### Feat

- **species**: added spin parsing from spin parity string

## v1.1.0 (2025-06-28)

### Feat

- **species**: added half life, spin parity, and decay modes to species database

## v1.0.9 (2025-06-26)

### Refactor

- **logs**: register symbol log info -> trace_l3

## v1.0.8 (2025-06-25)

### Feat

- **Composition**: added getMolarAbundance method

## v1.0.7 (2025-06-22)

### Fix

- **subprojecst**: version bump on liblogging and libconfig

## v1.0.6 (2025-06-22)

### Fix

- **headers**: moved all headers to fourdst/

## v1.0.5 (2025-06-21)

### Fix

- **header**: version

## v1.0.4 (2025-06-21)

### Fix

- **header**: fixed

## v1.0.3 (2025-06-21)

### Feat

- **header**: moved

## v1.0.2 (2025-06-21)

### Feat

- **header**: updated header location

## v1.0.1 (2025-06-21)

### Fix

- **subprojects**: removed quill and yaml-cpp meson build artifacts

## v1.0.0 (2025-06-21)

### Feat

- **network**: major progress on network finalization and matrix creation
- **assets/static**: moved data type logic to dedicated headers
- **composition**: added contains method
- **reaclib**: working on building efficient reaclib tooling for general nuclear network
- **reaclib**: working on building general, reaclib, based nuclear network
- **eos**: EOS now uses composition module
- **pythonInterface/eos**: fast forward
- **python-composition**: added composition module interface
- **pybind11**: added infra to compile with pybind11
- **composition**: added mix method to combine compositions. Also overloaded the + operator to mix with an assumed fraction of 50/50
- **composition**: added ability to change composition modes
- **composition**: added numberFrac methods and subset method
- **atomicSpecies.h**: regenerated with copy constructor
- **composition**: added composition class
- **composition**: added composition module stub
- **atomic-weights**: added AME2020 atomic masses
- **mfem**: added mfem into source tree along with patch based build system

### Fix

- **atomicSpecies.h-->-species.h**: added species.h includes for spesific species where needed
- **composition**: updated includes to include new assets/static/atomic header
- **src**: updated to compile on gcc and clang
- **eos**: fixed calculation of mean atomic number
- **eos-bindings**: minor bug fixes to bring eos bindings up to main
- **composition**: removed old py structure

### Refactor

- **liblogging**: changed SERiF to use liblogging
- **network**: updated network and network::approx8 to use composition module
- **serif**: updated tests to reflect new serif namespaces
- **serif**: refactored entire codebase into serif and sub namespaces
