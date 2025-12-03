logging_p = subproject('liblogging', default_options: [
  'pkg_config=' + get_option('pkg_config').to_string(),
  'build_tests=' + get_option('build_tests').to_string(),
  'build_examples=' + get_option('build_examples').to_string()
])

logging_dep  = logging_p.get_variable('logging_dep')
quill_dep = logging_p.get_variable('quill_dep')

log_dep = [logging_dep, quill_dep]
