config_p = subproject('libconfig', default_options: [
  'pkg_config=' + get_option('pkg_config').to_string(),
  'build_tests=' + get_option('build_tests').to_string(),
  'build_examples=' + get_option('build_examples').to_string()
  ])
config_dep  = config_p.get_variable('config_dep')
