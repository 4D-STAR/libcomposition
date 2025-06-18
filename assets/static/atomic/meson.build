required_headers = [
    'atomicSpecies.h',
]

foreach h : required_headers
    if not cpp.has_header(h, include_directories: include_directories('include'))
        error('SERiF requires the header file ' + h + ' to be present in the assets/static/atomic/include directory.')
    endif
endforeach

species_weight_dep = declare_dependency(
    include_directories: include_directories('include'),
)

message('✅ SERiF species_weight dependency declared')