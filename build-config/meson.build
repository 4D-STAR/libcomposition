cmake = import('cmake')

subdir('mfem')
subdir('yaml-cpp')
subdir('quill')
subdir('boost')
subdir('opatIO')
subdir('mpi')
subdir('hypre')
subdir('pybind')
subdir('cppad')

# Set the config file error handling options
configErr = get_option('config_error_handling')

# build up any -D flags we need
commonCppArgs = []
if configErr == 'warn'
    commonCppArgs += ['-DCONFIG_WARN']
elif configErr == 'harsh'
    commonCppArgs += ['-DCONFIG_HARSH']
endif

add_project_arguments(commonCppArgs, language: 'cpp')
