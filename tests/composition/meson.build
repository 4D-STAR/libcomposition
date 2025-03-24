# Test files for const
test_sources = [
    'compositionTest.cpp',
]

foreach test_file : test_sources
  exe_name = test_file.split('.')[0]
  message('Building test: ' + exe_name)

  # Create an executable target for each test
  test_exe = executable(
      exe_name,
      test_file,
      dependencies: [gtest_dep, species_weight_dep, gtest_main],
      install_rpath: '@loader_path/../../src'  # Ensure runtime library path resolves correctly
  )

  # Add the executable as a test
  test(
    exe_name,
    test_exe,
    env: ['MESON_SOURCE_ROOT=' + meson.project_source_root(), 'MESON_BUILD_ROOT=' + meson.project_build_root()])
endforeach
