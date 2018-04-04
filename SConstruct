# This is a comment


import os
import subprocess



AddOption('--cxx',
	dest='cxx',
	type='string',
	nargs=1,
	action='store',
	help='cxx compiler')


AddOption('--gprof',
	action="store_true",
	dest='gprof',
	default=False)

AddOption('--use_ceres',
        action="store_true",
        dest='use_ceres',
        default=False)

os.system('perl src/setup.pl')

mpi_mode = False

if GetOption('cxx') == "mpic++" or GetOption('cxx') == "mpicc":
	mpi_mode = True

include = [Dir('/opt/local/include/'), '#/src/external/include/']
lib_path = Dir('/opt/local/lib/');

libs = ['boost_system-mt', 'boost_filesystem-mt', 'boost_iostreams-mt' , 'boost_thread-mt', 'pthread', 'boost_serialization-mt', 'z']
if mpi_mode == True:
	libs = libs + ['boost_mpi-mt']


#svn_rev = os.popen('svnversion -n .').readline().rstrip()
#svn_rev = '\\"' + svn_rev + '\\"'
#hg_rev = os.popen('hg id').readline().rstrip()
#hg_rev = '\\"' + hg_rev + '\\"'
#pd2_ver = os.popen('cat src/Version.txt').readline().rstrip()
#pd2_ver = '\\"' + pd2_ver + '\\"' 
#comp_date = os.popen('date').readline().rstrip()
#comp_date = '\\"' + comp_date + '\\"'


cppdefines = ['USING_SCONS_BUILDER', ('HAVE_INLINE', 1), 'NDEBUG', ('BOOST_FILESYSTEM_VERSION', 3) ] #[('SVN_REV',svn_rev), ('HG_REV',hg_rev), ('_PRODART_VERSION_', pd2_ver), ('_COMPILE_DATE_', comp_date)]
if mpi_mode == True:
	cppdefines = cppdefines + ['PD2_MPI_MODE']

gcc_cppflags = ['-ffast-math', '-msse3', '-funroll-loops', '-pipe','-g', '-Wall', '-fmessage-length=0', '-std=c++11' ]

# optimisation level normally -O2
gcc_cppflags = gcc_cppflags + ['-O3']

# link time optimisation
#gcc_cppflags = gcc_cppflags + ['-flto']


#gcc_cppflags = gcc_cppflags + ['-rdynamic']

# FOR profiling
if GetOption('gprof') == True:
	gcc_cppflags = gcc_cppflags + ['-pg']
	print "INFO: compiling with gprof flags"

if GetOption('use_ceres') == True:
        cppdefines = cppdefines + ['PD2_USE_CERES']
	libs = libs + ['ceres'] + ['glog', 'cxsparse', 'cholmod', 'lapack', 'blas']
#	libs = libs + ['SuiteSparse']			# this doesn't seem to exist on RHEL 7 or Oracle linux 7 
#	include = include + [Dir('/opt/local/include/eigen3/')]
        print "INFO: compiling and linking to use Ceres-solver"
else:
	libs = libs + ['gsl', 'gslcblas']
	cppdefines = cppdefines + [('GSL_RANGE_CHECK', 0)]

#pose_src = Glob('src/pose//*.cpp')

bbb_src = Glob('src/backbonebuilder/*.cpp')

movers_src = Glob('src/movers/*.cpp') +  Glob('src/movers/bb_movers/*.cpp') + Glob('src/movers/ca_movers/*.cpp') +  Glob('src/movers/fa_movers/*.cpp')

pose_src = Glob('src/pose/*.cpp')

pose_meta_src = Glob('src/pose_meta/*.cpp') +  Glob('src/pose_meta/frag_classify/*.cpp')

pose_utils_src = Glob('src/pose_utils/*.cpp')

potentials_src = Glob('src/potentials/*.cpp') +  Glob('src/potentials/bb_potentials/*.cpp') +  Glob('src/potentials/ca_potentials/*.cpp')

pd_env_src = Glob('src/prodart_env/*.cpp')

protocols_src = Glob('src/protocols/*.cpp')
mpi_protocols_src = Glob('src/protocols/mpi/*.cpp')
rotamers_src = Glob('src/rotamers/*.cpp') + Glob('src/rotamers/sidechain_builders/*.cpp')

sim_src = Glob('src/simulation/*.cpp') +  Glob('src/simulation/bb_mc_protocols/*.cpp') + Glob('src/simulation/ca_mc_protocols/*.cpp')
mpi_sim_src = Glob('src/simulation/mpi/*.cpp') +  Glob('src/simulation/bb_mc_protocols/mpi/*.cpp') + Glob('src/simulation/ca_mc_protocols/mpi/*.cpp')
sticks_src = Glob('src/sticks/*.cpp')

tclap_src = Glob('src/tclap/*.cpp')

utils_src = Glob('src/utils/*.cpp')


common_sources = bbb_src + movers_src + pose_src + pose_meta_src + pose_utils_src + potentials_src + pd_env_src + protocols_src + rotamers_src + sim_src + sticks_src + tclap_src + utils_src
mpi_common_sources = mpi_sim_src + mpi_protocols_src

#print common_sources

gcc_link_flags = []
#gcc_link_flags = gcc_link_flags + ['-static']
gcc_link_flags = gcc_link_flags + ['-g']
gcc_link_flags = gcc_link_flags + ['-rdynamic']
#gcc_link_flags = gcc_link_flags + ['-flto']
#gcc_link_flags = gcc_link_flags + ['-Wl,-dead_strip']


if GetOption('gprof') == True:
	gcc_link_flags = gcc_link_flags + ['-pg']
	print "INFO: linking with gprof flags"


# for profiling:
#gcc_link_flags = gcc_link_flags + ['-pg']

gcc_env = Environment(ENV = os.environ, CC = 'gcc', CXX = 'g++', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags, LINKFLAGS=gcc_link_flags )   # Create an environmnet
gcc_env.Append(CPPPATH=["#/src"] + [include])

icc_env = Environment(ENV = os.environ, CC = 'icc', CXX = 'icpc', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags + ['-ip'], LINKFLAGS=gcc_link_flags + ['-ip'] )
icc_env.Append(CPPPATH=["#/src"] + [include])

#llvm_gcc_env = Environment(ENV = os.environ, CC = 'llvm-gcc', CXX = 'llvm-g++', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags, LINKFLAGS=gcc_link_flags )   # Create an environmnet
#llvm_gcc_env.Append(CPPPATH=["#/src"] + [include])

clang_env = Environment(ENV = os.environ, CC = 'clang', CXX = 'clang++', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags, LINKFLAGS=gcc_link_flags )   # Create an environmnet
clang_env.Append(CPPPATH=["#/src"] + [include])

#mpi_env = Environment(ENV = os.environ, CC = 'mpicc', CXX = 'mpic++', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags, LINKFLAGS=gcc_link_flags ) 
#mpi_env.Append(CPPPATH=["#/src"] + [include])


mpi_mode = False


if GetOption('cxx') == None:
	print "--cxx option not set: ", GetOption('cxx')
	env = gcc_env
	print "setting env to gcc_env"
elif GetOption('cxx') == "g++" or GetOption('cxx') == "gcc":
        print "--cxx option set: ", GetOption('cxx')
        env = gcc_env
        print "setting env to gcc_env"
elif GetOption('cxx') == "llvm-g++" or GetOption('cxx') == "llvm-gcc":
        print "--cxx option set: ", GetOption('cxx')
	env = llvm_gcc_env
	print "setting env to llvm_gcc_env"
elif GetOption('cxx') == "clang++" or GetOption('cxx') == "clang":
	print "--cxx option set: ", GetOption('cxx')
	env = clang_env
	print "setting env to clang_env"
elif GetOption('cxx') == "icc" or GetOption('cxx') == "icpc":
        print "--cxx option set: ", GetOption('cxx')
        env = icc_env
        print "setting env to icc_env"
elif GetOption('cxx') == "mpic++" or GetOption('cxx') == "mpicc":
        print "--cxx option set: ", GetOption('cxx')
        env = mpi_env
        print "setting env to mpi_env"
	mpi_mode = True
	#libs = libs + ['boost_mpi-mt']
	#mpi_env = Environment(ENV = os.environ, CC = 'mpicc', CXX = 'mpic++', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags, LINKFLAGS=gcc_link_flags )
	#mpi_env.Append(CPPPATH=["#/src"] + [include])
else:
        print "--cxx option set: ", GetOption('cxx')
	print "ERROR: C++ compiler option not recognised"



print "CC is:", env['CC']
print "CXX is:", env['CXX']
print "LINK is:", env['LINK']
print "TOOLS is:", env['TOOLS']


#print env.Dump()

try:
	os_cppath = os.environ['INCLUDE'].split(':')
	env.Append(CPPPATH=os_cppath)
except KeyError:
	print "INFO: no OS INCLUDE path found"	

print "CPPPATH =", env['CPPPATH']

try:
	os_libpath = os.environ['LD_LIBRARY_PATH'].split(':')
	env.Append(LIBPATH=os_libpath)
except KeyError:
	print "INFO: no OS LIBPATH path found"

print "LIBPATH =", env['LIBPATH']

if mpi_mode == False:
	env.Program(target = "bin/pd2_loop_model", source = ["src/pd2_loop_model.cpp"] + common_sources )
	env.Program(target = "bin/pd2_ca2main", source = ["src/pd2_ca2main.cpp"] + common_sources )
	env.Program(target = "bin/pd2_minimise", source = ["src/pd2_minimise.cpp"] + common_sources )
	env.Program(target = "bin/pdb_utils", source = ["src/pdb_utils.cpp"] + common_sources )
	env.Program(target = "bin/pd2_train", source = ["src/pd2_train.cpp"] + common_sources )
	#env.Program(target = "bin/", source = ["src/"] + common_sources )

else:
	print "INFO: mpi_mode"





