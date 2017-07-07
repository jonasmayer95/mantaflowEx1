file(REMOVE_RECURSE
  "pp/source/general.cpp"
  "pp/source/fluidsolver.cpp"
  "pp/source/conjugategrad.cpp"
  "pp/source/multigrid.cpp"
  "pp/source/grid.cpp"
  "pp/source/grid4d.cpp"
  "pp/source/levelset.cpp"
  "pp/source/fastmarch.cpp"
  "pp/source/shapes.cpp"
  "pp/source/mesh.cpp"
  "pp/source/particle.cpp"
  "pp/source/movingobs.cpp"
  "pp/source/fileio.cpp"
  "pp/source/noisefield.cpp"
  "pp/source/kernel.cpp"
  "pp/source/vortexsheet.cpp"
  "pp/source/vortexpart.cpp"
  "pp/source/turbulencepart.cpp"
  "pp/source/timing.cpp"
  "pp/source/edgecollapse.cpp"
  "pp/source/plugin/advection.cpp"
  "pp/source/plugin/extforces.cpp"
  "pp/source/plugin/flip.cpp"
  "pp/source/plugin/fire.cpp"
  "pp/source/plugin/fluidguiding.cpp"
  "pp/source/plugin/kepsilon.cpp"
  "pp/source/plugin/initplugins.cpp"
  "pp/source/plugin/meshplugins.cpp"
  "pp/source/plugin/pressure.cpp"
  "pp/source/plugin/surfaceturbulence.cpp"
  "pp/source/plugin/vortexplugins.cpp"
  "pp/source/plugin/waveletturbulence.cpp"
  "pp/source/plugin/waves.cpp"
  "pp/source/python/defines.py"
  "pp/source/python/defines.py.reg"
  "pp/source/test.cpp"
  "pp/source/ex1.cpp"
  "pp/source/ex2.cpp"
  "pp/source/ex3.cpp"
  "pp/source/ex4.cpp"
  "pp/source/gui/customctrl.cpp"
  "pp/source/gui/mainwindow.cpp"
  "pp/source/gui/glwidget.cpp"
  "pp/source/gui/painter.cpp"
  "pp/source/gui/meshpainter.cpp"
  "pp/source/gui/particlepainter.cpp"
  "pp/source/gui/qtmain.cpp"
  "pp/source/general.h"
  "pp/source/general.h.reg"
  "pp/source/commonkernels.h"
  "pp/source/commonkernels.h.reg"
  "pp/source/conjugategrad.h"
  "pp/source/conjugategrad.h.reg"
  "pp/source/multigrid.h"
  "pp/source/multigrid.h.reg"
  "pp/source/fastmarch.h"
  "pp/source/fastmarch.h.reg"
  "pp/source/fluidsolver.h"
  "pp/source/fluidsolver.h.reg"
  "pp/source/grid.h"
  "pp/source/grid.h.reg"
  "pp/source/grid4d.h"
  "pp/source/grid4d.h.reg"
  "pp/source/mesh.h"
  "pp/source/mesh.h.reg"
  "pp/source/particle.h"
  "pp/source/particle.h.reg"
  "pp/source/levelset.h"
  "pp/source/levelset.h.reg"
  "pp/source/shapes.h"
  "pp/source/shapes.h.reg"
  "pp/source/noisefield.h"
  "pp/source/noisefield.h.reg"
  "pp/source/vortexsheet.h"
  "pp/source/vortexsheet.h.reg"
  "pp/source/kernel.h"
  "pp/source/kernel.h.reg"
  "pp/source/timing.h"
  "pp/source/timing.h.reg"
  "pp/source/movingobs.h"
  "pp/source/movingobs.h.reg"
  "pp/source/fileio.h"
  "pp/source/fileio.h.reg"
  "pp/source/edgecollapse.h"
  "pp/source/edgecollapse.h.reg"
  "pp/source/vortexpart.h"
  "pp/source/vortexpart.h.reg"
  "pp/source/turbulencepart.h"
  "pp/source/turbulencepart.h.reg"
  "pp/source/pcgsolver.h"
  "pp/source/pcgsolver.h.reg"
  "pp/source/gui/mainwindow.h"
  "pp/source/gui/mainwindow.h.reg"
  "pp/source/gui/glwidget.h"
  "pp/source/gui/glwidget.h.reg"
  "pp/source/gui/painter.h"
  "pp/source/gui/painter.h.reg"
  "pp/source/gui/meshpainter.h"
  "pp/source/gui/meshpainter.h.reg"
  "pp/source/gui/qtmain.h"
  "pp/source/gui/qtmain.h.reg"
  "pp/source/gui/customctrl.h"
  "pp/source/gui/customctrl.h.reg"
  "pp/source/gui/particlepainter.h"
  "pp/source/gui/particlepainter.h.reg"
  "pp/source/python/defines.py.reg.cpp"
  "pp/source/general.h.reg.cpp"
  "pp/source/commonkernels.h.reg.cpp"
  "pp/source/conjugategrad.h.reg.cpp"
  "pp/source/multigrid.h.reg.cpp"
  "pp/source/fastmarch.h.reg.cpp"
  "pp/source/fluidsolver.h.reg.cpp"
  "pp/source/grid.h.reg.cpp"
  "pp/source/grid4d.h.reg.cpp"
  "pp/source/mesh.h.reg.cpp"
  "pp/source/particle.h.reg.cpp"
  "pp/source/levelset.h.reg.cpp"
  "pp/source/shapes.h.reg.cpp"
  "pp/source/noisefield.h.reg.cpp"
  "pp/source/vortexsheet.h.reg.cpp"
  "pp/source/kernel.h.reg.cpp"
  "pp/source/timing.h.reg.cpp"
  "pp/source/movingobs.h.reg.cpp"
  "pp/source/fileio.h.reg.cpp"
  "pp/source/edgecollapse.h.reg.cpp"
  "pp/source/vortexpart.h.reg.cpp"
  "pp/source/turbulencepart.h.reg.cpp"
  "pp/source/pcgsolver.h.reg.cpp"
  "pp/source/gui/mainwindow.h.reg.cpp"
  "pp/source/gui/glwidget.h.reg.cpp"
  "pp/source/gui/painter.h.reg.cpp"
  "pp/source/gui/meshpainter.h.reg.cpp"
  "pp/source/gui/qtmain.h.reg.cpp"
  "pp/source/gui/customctrl.h.reg.cpp"
  "pp/source/gui/particlepainter.h.reg.cpp"
  "pp/source/registration.cpp"
  "pp/source/gui/moc_mainwindow.cpp"
  "pp/source/gui/moc_glwidget.cpp"
  "pp/source/gui/moc_painter.cpp"
  "pp/source/gui/moc_meshpainter.cpp"
  "pp/source/gui/moc_qtmain.cpp"
  "pp/source/gui/moc_customctrl.cpp"
  "pp/source/gui/moc_particlepainter.cpp"
  "qrc_res.cpp"
  "pp/source/gitinfo.h"
  "CMakeFiles/manta.dir/pwrapper/pymain.cpp.o"
  "CMakeFiles/manta.dir/pwrapper/pclass.cpp.o"
  "CMakeFiles/manta.dir/pwrapper/pvec3.cpp.o"
  "CMakeFiles/manta.dir/pwrapper/pconvert.cpp.o"
  "CMakeFiles/manta.dir/pwrapper/registry.cpp.o"
  "CMakeFiles/manta.dir/util/vectorbase.cpp.o"
  "CMakeFiles/manta.dir/util/vector4d.cpp.o"
  "CMakeFiles/manta.dir/util/simpleimage.cpp.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/adler32.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/compress.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/crc32.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/deflate.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/gzclose.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/gzlib.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/gzread.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/gzwrite.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/inflate.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/infback.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/inftrees.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/inffast.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/trees.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/uncompr.c.o"
  "CMakeFiles/manta.dir/dependencies/zlib-1.2.8/zutil.c.o"
  "CMakeFiles/manta.dir/pp/source/general.cpp.o"
  "CMakeFiles/manta.dir/pp/source/fluidsolver.cpp.o"
  "CMakeFiles/manta.dir/pp/source/conjugategrad.cpp.o"
  "CMakeFiles/manta.dir/pp/source/multigrid.cpp.o"
  "CMakeFiles/manta.dir/pp/source/grid.cpp.o"
  "CMakeFiles/manta.dir/pp/source/grid4d.cpp.o"
  "CMakeFiles/manta.dir/pp/source/levelset.cpp.o"
  "CMakeFiles/manta.dir/pp/source/fastmarch.cpp.o"
  "CMakeFiles/manta.dir/pp/source/shapes.cpp.o"
  "CMakeFiles/manta.dir/pp/source/mesh.cpp.o"
  "CMakeFiles/manta.dir/pp/source/particle.cpp.o"
  "CMakeFiles/manta.dir/pp/source/movingobs.cpp.o"
  "CMakeFiles/manta.dir/pp/source/fileio.cpp.o"
  "CMakeFiles/manta.dir/pp/source/noisefield.cpp.o"
  "CMakeFiles/manta.dir/pp/source/kernel.cpp.o"
  "CMakeFiles/manta.dir/pp/source/vortexsheet.cpp.o"
  "CMakeFiles/manta.dir/pp/source/vortexpart.cpp.o"
  "CMakeFiles/manta.dir/pp/source/turbulencepart.cpp.o"
  "CMakeFiles/manta.dir/pp/source/timing.cpp.o"
  "CMakeFiles/manta.dir/pp/source/edgecollapse.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/advection.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/extforces.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/flip.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/fire.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/fluidguiding.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/kepsilon.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/initplugins.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/meshplugins.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/pressure.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/surfaceturbulence.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/vortexplugins.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/waveletturbulence.cpp.o"
  "CMakeFiles/manta.dir/pp/source/plugin/waves.cpp.o"
  "CMakeFiles/manta.dir/pp/source/test.cpp.o"
  "CMakeFiles/manta.dir/pp/source/ex1.cpp.o"
  "CMakeFiles/manta.dir/pp/source/ex2.cpp.o"
  "CMakeFiles/manta.dir/pp/source/ex3.cpp.o"
  "CMakeFiles/manta.dir/pp/source/ex4.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/customctrl.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/mainwindow.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/glwidget.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/painter.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/meshpainter.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/particlepainter.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/qtmain.cpp.o"
  "CMakeFiles/manta.dir/pp/source/python/defines.py.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/general.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/commonkernels.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/conjugategrad.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/multigrid.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/fastmarch.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/fluidsolver.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/grid.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/grid4d.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/mesh.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/particle.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/levelset.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/shapes.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/noisefield.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/vortexsheet.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/kernel.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/timing.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/movingobs.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/fileio.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/edgecollapse.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/vortexpart.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/turbulencepart.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/pcgsolver.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/mainwindow.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/glwidget.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/painter.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/meshpainter.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/qtmain.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/customctrl.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/particlepainter.h.reg.cpp.o"
  "CMakeFiles/manta.dir/pp/source/registration.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/moc_mainwindow.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/moc_glwidget.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/moc_painter.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/moc_meshpainter.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/moc_qtmain.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/moc_customctrl.cpp.o"
  "CMakeFiles/manta.dir/pp/source/gui/moc_particlepainter.cpp.o"
  "CMakeFiles/manta.dir/qrc_res.cpp.o"
  "manta.pdb"
  "manta"
)

# Per-language clean rules from dependency scanning.
foreach(lang C CXX)
  include(CMakeFiles/manta.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
