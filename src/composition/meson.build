required_headers = [
    'fourdst/composition/atomicSpecies.h',
    'fourdst/composition/species.h',
    'fourdst/composition/composition.h'
]

foreach h : required_headers
    if not cpp.has_header(h, include_directories: include_directories('include'))
        error('libcomposition requires the header file ' + h + ' to be present in the fourdst/composition directory.')
    endif
endforeach

species_weight_dep = declare_dependency(
    include_directories: include_directories('include'),
)

message('✅ libcomposition species_weight dependency declared')

composition_sources = files(
  'lib/composition.cpp',
)


dependencies = [
    species_weight_dep,
    const_dep,
    config_dep,
    log_dep
]

# Define the libcomposition library so it can be linked against by other parts of the build system
libcomposition = library('composition',
    composition_sources, 
    cpp_args: ['-fvisibility=default'],
    dependencies: dependencies,
    install : true)

composition_dep = declare_dependency(
    link_with: libcomposition,
    dependencies: dependencies,
    sources: composition_sources,
)

# Make headers accessible
composition_headers = files(
  'include/fourdst/composition/composition.h',
  'include/fourdst/composition/atomicSpecies.h',
  'include/fourdst/composition/species.h',
  'include/fourdst/composition/elements.h'
)
install_headers(composition_headers, subdir : 'fourdst/fourdst/composition')

composition_exception_headers = files(
    'include/fourdst/composition/exceptions/exceptions_composition.h',
)
install_headers(composition_exception_headers, subdir : 'fourdst/fourdst/composition/exceptions')
