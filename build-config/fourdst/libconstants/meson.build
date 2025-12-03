const_p = subproject(
  'libconstants',
  default_options: [
    'pkg_config=' + get_option('pkg_config').to_string(),
    'build_tests=' + get_option('build_tests').to_string(),
    'build_examples=' + get_option('build_examples').to_string()
]
)
const_dep  = const_p.get_variable('const_dep')
