# Build list of app dependencies
app_deps = [moist_dep]
app_deps += fclap_dep

moist_exe = executable(
   meson.project_name(),
   sources: files(
     'main.f90',
     'cli.f90',
     'driver.f90',
   ),
   dependencies: app_deps,
   install: install,
   link_language: 'fortran',
 )

app_tester = find_program(files('tester.py'))
