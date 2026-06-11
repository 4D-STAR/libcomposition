# ***********************************************************************
#
#   Copyright (C) 2025 -- The 4D-STAR Collaboration
#   File Author: Emily Boudreaux
#   Last Modified: June 21, 2025
#
#   libcomposition is free software; you can use it and/or modify
#   it under the terms and restrictions the GNU General Library Public
#   License version 3 (GPLv3) as published by the Free Software Foundation.
#
#   libcomposition is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#   See the GNU Library General Public License for more details.
#
#   You should have received a copy of the GNU Library General Public License
#   along with this software; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# *********************************************************************** #
project('libcomposition', 'cpp', version: 'v2.4.6', default_options: ['cpp_std=c++23'], meson_version: '>=1.5.0')

# Add default visibility for all C++ targets
add_project_arguments('-fvisibility=default', language: 'cpp')
# Disable shadow warnings
add_project_arguments('-Wno-shadow', language: 'cpp')

if get_option('build_python')
    local_py_install = import('python').find_installation('python3', pure: false)
    fourdst_pkg = local_py_install.get_install_dir() / 'fourdst'

    composition_libdir             = fourdst_pkg / 'lib'
    composition_header_install_dir = fourdst_pkg / 'include' / 'fourdst' / 'composition'

    vendor_lib_dir             = fourdst_pkg / 'lib' / 'vendor'
    vendor_header_dir          = fourdst_pkg / 'include' / 'fourdst' / 'vendor'
    message('Installing libcomposition in python mode. Install locations are...')
else
    composition_libdir             = get_option('libdir')
    composition_header_install_dir = get_option('includedir') / 'fourdst' / 'composition'
    message('Installing libcomposition in non python (normal) mode. Install locations are...')
endif
message('    Libraries: ' + composition_libdir)
message('    Headers:   ' + composition_header_install_dir)

cpp = meson.get_compiler('cpp')
subdir('build-config')


subdir('src')

if get_option('build_tests')
  subdir('tests')
endif

if get_option('build_examples')
  subdir('examples')
endif

if get_option('build_benchmarks')
    subdir('benchmarks')
endif

if get_option('pkg_config')
    message('Generating pkg-config file for libcomposition...')
    pkg = import('pkgconfig')
    pkg.generate(
        name: 'libcomposition',
        description: 'Composition module for SERiF and related projects',
        version: meson.project_version(),
        libraries: [libcomposition, '-Wl,-rpath,${libdir}'],
        filebase: 'fourdst_composition',
        install_dir: join_paths(get_option('libdir'), 'pkgconfig')
    )
endif