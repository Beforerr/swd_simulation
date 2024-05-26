vpic := home_dir + "/src/vpic-kokkos/build/bin/vpic"
project := "UCLA0040"
queue := "main"

compile:
  cd {{invocation_directory()}}; {{vpic}} *.cxx

q-sub:
  cd {{invocation_directory()}}; qsub -A {{project}} -q {{queue}}  *.pbs

q-subf file:
  cd {{invocation_directory()}}; qsub -A {{project}} -q {{queue}}  {{file}}

q-info:
  qhist -u $USER
  -qstat -u $USER

vpic:
  micromamba env create vtk pyvista pyqt VisualPIC --name vpic
  micromamba run -n vpic vpic -h
  pipx install VisualPIC --preinstall pyqt5 --preinstall vtk --preinstall pyvista
  pipx inject visualpic pyqt5 vtk pyvista
  vpic dim_2_beta_0.25_theta_60/diags/diag1/ -Bx -By -Bz
  vpic dim_2_beta_0.25_theta_60/diags/diag1/ -Jx -Jy -Jz

  micromamba run -n vpic vpic3d dim_3_beta_0.25_theta_60/diags/diag1/ -Bx -By -Bz
  micromamba run -n vpic vpic3d dim_3_beta_0.25_theta_60/diags/diag1/ -Jx
  micromamba run -n vpic vpic3d dim_3_beta_0.25_theta_60/diags/diag1/ -Jz