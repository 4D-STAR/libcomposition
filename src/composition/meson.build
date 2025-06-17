# Define the library
composition_sources = files(
  'private/composition.cpp',
)

composition_headers = files(
    'public/composition.h'
)

dependencies = [
    probe_dep,
    quill_dep,
    species_weight_dep
]

# Define the libcomposition library so it can be linked against by other parts of the build system
libcomposition = library('composition',
    composition_sources, 
    include_directories: include_directories('public'),
    cpp_args: ['-fvisibility=default'],
    dependencies: dependencies,
    install : true)

composition_dep = declare_dependency(
    include_directories: include_directories('public'),
    link_with: libcomposition,
    dependencies: dependencies,
    sources: composition_sources,
)

# Make headers accessible
install_headers(composition_headers, subdir : 'SERiF/composition')