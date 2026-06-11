required_headers = [
    'fourdst/atomic/atomicSpecies.h',
    'fourdst/atomic/species.h',
    'fourdst/composition/composition.h',
    'fourdst/composition/utils/utils.h',
    'fourdst/composition/composition_abstract.h',
    'fourdst/composition/exceptions/exceptions_composition.h',
    'fourdst/composition/io/standard_compositions.h',
    'fourdst/composition/io/StandardMetalFractionsBinary.h'
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
  'lib/utils.cpp',
  'lib/decorators/composition_masked.cpp',
  'lib/io/standard_compositions.cpp'
)


dependencies = [
    species_weight_dep,
    const_dep,
    config_dep,
    log_dep,
    xxhash_dep
]

fourdst_component_libdir = get_option('lib_install_dir') != '' ? get_option('lib_install_dir') : get_option('libdir')

libcomposition = library('composition',
                         composition_sources,
                         cpp_args: ['-fvisibility=default'],
                         dependencies: dependencies,
                         install: true,
                         install_dir: fourdst_component_libdir)


composition_dep = declare_dependency(
    link_with: libcomposition,
    dependencies: dependencies,
    sources: composition_sources,
)

# Make headers accessible
composition_headers = files(
  'include/fourdst/composition/composition.h',
  'include/fourdst/composition/composition_abstract.h',
)
install_headers(composition_headers, subdir : 'fourdst/composition')

composition_headers_utils = files(
    'include/fourdst/composition/utils/utils.h',
    'include/fourdst/composition/utils/composition_hash.h'
)
install_headers(composition_headers_utils, subdir : 'fourdst/composition/utils')

composition_headers_io = files(
    'include/fourdst/composition/io/standard_compositions.h',
    'include/fourdst/composition/io/StandardMetalFractionsBinary.h'
)

install_headers(composition_headers_io, subdir : 'fourdst/composition/io')

composition_headers_decorators = files(
    'include/fourdst/composition/decorators/composition_masked.h',
    'include/fourdst/composition/decorators/composition_decorator_abstract.h',
)

install_headers(composition_headers_decorators, subdir : 'fourdst/composition/decorators')

composition_headers_atomic = files(
    'include/fourdst/atomic/atomicSpecies.h',
    'include/fourdst/atomic/elements.h',
    'include/fourdst/atomic/species.h',
    
)
install_headers(composition_headers_atomic, subdir : 'fourdst/atomic')

composition_exception_headers = files(
    'include/fourdst/composition/exceptions/exceptions_composition.h',
)
install_headers(composition_exception_headers, subdir : 'fourdst/composition/exceptions')

composition_iterator_headers = files(
    'include/fourdst/composition/iterators/composition_abstract_iterator.h',
)

install_headers(composition_iterator_headers, subdir : 'fourdst/composition/iterators')

v = meson.project_version()

conf_data = configuration_data()
conf_data.set_quoted('PROJECT_VERSION', v)
conf_data.set_quoted('PROJECT_NAME', meson.project_name())

configure_file(
    input : 'include/fourdst/config.h.in',
    output : 'config.h',
    configuration : conf_data
)
