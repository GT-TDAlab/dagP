#-*-python-*-

import string
import os
import sys


## Make it colorful ......

colors = {
    'white':    "\033[1;37m",
    'yellow':   "\033[1;33m",
    'green':    "\033[1;32m",
    'blue':     "\033[1;34m",
    'cyan':     "\033[1;36m",
    'red':      "\033[1;31m",
    'purple':   "\033[95m",
    'magenta':  "\033[1;35m",
    'black':      "\033[1;30m",
    'darkwhite':  "\033[0;37m",
    'darkyellow': "\033[0;33m",
    'darkgreen':  "\033[0;32m",
    'darkblue':   "\033[0;34m",
    'darkcyan':   "\033[0;36m",
    'darkred':    "\033[0;31m",
    'darkmagenta':"\033[0;35m",
    'darkblack':  "\033[0;30m",
    'end':        "\033[0;0m"
}


#If the output is not a terminal, remove the colors
if not sys.stdout.isatty():
   for key, value in colors.items():
      colors[key] = ''



def DumpEnv( env, key = None, header = None, footer = None ):
  import pprint
  pp = pprint.PrettyPrinter( indent = 2 )
  if key:
     dict = env.Dictionary( key )
  else:
     dict = env.Dictionary()
     if header:
        print(header)
     pp.pprint( dict )
     if footer:
        print(footer)


vars = Variables('config.py')
vars.AddVariables(
  ('metis', "set to 1 to build with metis", 0),
  ('scotch', "set to 1 to build with scotch", 0),
  BoolVariable('debug', "set to 1 to build in debug mode", 0),
  ('debug_flags', "additional debug flags", None),
  ('build_threads', "additional debug flags", 1),
  BoolVariable('m64', "set to 1 to build 64-bit", 1),
  ('CC', "set C compiler", 'cc'),
  ('CXX', 'Forces C++ compiler', None),
  ('defines', 'List any additional preprocessor defines here', ''),
  ('CCFLAGS', 'Forces C compiler flags', None),
  ('CXXFLAGS', 'Forces C++ compiler flags', None),
  ('LINK', 'Forces linker ', None),
  ('LINKFLAGS', 'Forces linker flags', None),
  ('extincludes', 'List any additional include paths for extensions here', None),
  ('extlibs', 'List any additional link libraries for extensions here', None),
  ('extlibpath', 'List any additional link paths for extensions here', None),
  ('extdefines', 'List any additional preprocessor defines for extensions here', None),
)

envC = Environment(variables=vars, ENV = os.environ)

name = os.uname()[0]
arch = os.uname()[4]
print(name, arch)

force32v64 = 0
if (name == 'Darwin') or (name == 'Linux'):
   force32v64 = 1
   if (int(envC.get('m64', 1))==1):
      arch = 'x86_64'
   else:
      arch = 'i386'
platform = name + "-" + ''.join(arch.split())


print("%sBuilding for platform <%s" % (colors['red'], colors['darkred']), platform, "%s>%s" % (colors['red'], colors['end']))


compile_source_message = '%sCompiling %s==> %s$SOURCE%s' % \
   (colors['blue'], colors['purple'], colors['darkgreen'], colors['end'])

compile_shared_source_message = '%sCompiling shared %s==> %s$SOURCE%s' % \
   (colors['blue'], colors['purple'], colors['darkred'], colors['end'])

link_program_message = '%sLinking Program %s==> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['darkblack'], colors['end'])

link_library_message = '%sLinking Static Library %s==> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['darkmagenta'], colors['end'])

ranlib_library_message = '%sRanlib Library %s==> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['darkmagenta'], colors['end'])

link_shared_library_message = '%sLinking Shared Library %s==> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['darkmagenta'], colors['end'])

java_library_message = '%sCreating Java Archive %s==> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['darkgreen'], colors['end'])
tar_message = '%sCreating Tar Archive %s==> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['darkmagenta'], colors['end'])


include = "#/src/common:#/src/recBisection".split(':')
lib = "#/lib"
bin = "#/exe"

env = Environment(variables=vars, ENV = os.environ, TARFLAGS = '-c -z',
  CXXCOMSTR = compile_source_message + "\n--> $CXXCOM",
  CCCOMSTR = compile_source_message + "\n--> $CCCOM",
  SHCCCOMSTR = compile_shared_source_message+"\n--> $SHCCCOM",
  SHCXXCOMSTR = compile_shared_source_message+"\n--> $SHCXXCOM",
  ARCOMSTR = link_library_message+"\n--> $ARCOM",
  RANLIBCOMSTR = ranlib_library_message+"\n--> $RANLIBCOM",
  SHLINKCOMSTR = link_shared_library_message+"\n--> $SHLINKCOM",
  LINKCOMSTR = link_program_message+"\n--> $LINKCOM",
  JARCOMSTR = java_library_message+"\n--> $JARCOM",
  JAVACCOMSTR = compile_source_message+"\n--> $JAVACCOM",
  TARCOMSTR = tar_message+"\n--> $TARCOM",
  PLATFORM = platform,
  BINDIR = bin,
  INCDIR = include,
  LIBDIR = lib,
  CPPPATH = [include],
  LIBPATH = [lib]
)

env.Help(vars.GenerateHelpText(env))

SetOption('num_jobs', env.get("build_threads"))
print ("running with -j", GetOption('num_jobs'))

metis = int(env.get('metis', 0))
scotch = int(env.get('scotch', 0))
debug_mode = env.get('debug', 0)
debug_flags = env.get('debug_flags', 0)
defines_st = env.get('defines')

libs = []
defines = defines_st.split(' ')

if debug_mode:
    print("%sDebug%s build..." % (colors['green'],colors['end']))
    defines.append('_DEBUG')
    env['CCFLAGS'] = env['CCFLAGS'] + ' -g'
    if debug_flags:
        env['CCFLAGS'] = env['CCFLAGS'] + ' ' + debug_flags
else:
    print("%sRelease%s build..." % (colors['green'],colors['end']))
    env['CCFLAGS'] = env['CCFLAGS'] + ' -O3'

if metis == 1:
   defines.append('dagP_METIS')

if scotch == 1:
   defines.append('dagP_SCOTCH')

if debug_mode >= 2:
   defines.append('_DEBUG1')
if debug_mode >= 3:
   defines.append('_DEBUG2')
if debug_mode >= 4:
   defines.append('_DEBUG3')


if force32v64:
   m64 = int(env.get('m64', 1))
   if m64:
      env['CCFLAGS'] = env['CCFLAGS'] + ' -m64'
      env['CXXFLAGS'] = env['CXXFLAGS'] + ' -m64'
      env['LINKFLAGS'] = env['LINKFLAGS'] + ' -m64'
   else:
      env['CCFLAGS'] = env['CCFLAGS'] + ' -m32'
      env['CXXFLAGS'] = env['CXXFLAGS'] + ' -m32'
      env['LINKFLAGS'] = env['LINKFLAGS'] + ' -m32'

env['CXXFLAGS'] = env['CXXFLAGS'] + ' -std=c++11'


# add -fPIC for x86_64-bit
if os.uname()[4] == 'x86_64':
    env['CCFLAGS'] = env['CCFLAGS'] + ' -fPIC'
    env['CXXFLAGS'] = env['CXXFLAGS'] + ' -fPIC'

if name == 'Darwin':
   env['SHLINKFLAGS'] = '$LINKFLAGS -dynamic'
   env['SHLIBSUFFIX'] = '.dylib'
   env['CCFLAGS'] = env['CCFLAGS'] + ' -Wall -Wstrict-prototypes -Wno-unused-variable -Wno-unused-parameter  -Werror-implicit-function-declaration -W'

elif name == 'Linux':
   env['CCFLAGS'] = env['CCFLAGS'] + ' -Wall -Wstrict-prototypes -Wno-unused-variable  -Werror-implicit-function-declaration -W'
   env.Append(LINKFLAGS =['-static'])
   if debug_mode==0:
         env.Append(LINKFLAGS =['-s'])


env.Append(CPPDEFINES=defines)

env.Append(CPPPATH=['#/src'])
# env.Append(LIBPATH='$PLATFORM/lib')

dgraphlibs = ['dagp']
libs.extend(dgraphlibs);
env.Append(LIBS = libs)


dgraphlibname = 'libdagp.a'

env['dgraphlibname'] = dgraphlibname
env['dgraphlibs'] = dgraphlibs

Export("env")

###
### external libraries
###
extincludes = env.get('extincludes')
extlibs = env.get('extlibs')
extlibpath = env.get('extlibpath')
extdefines = env.get('extdefines')

if extincludes:
    env.Append(CPPPATH=extincludes.split(':'))
if extlibs:
    added = extlibs.split(' ')
    env.Append(LIBS = added)
if extlibpath:
    env.Append(LIBPATH=extlibpath.split(':'))
if extdefines:
    env.Append(CPPDEFINES=extdefines.split(' '))

env.Append(LIBS = ['m'])


dgraphlibsrcs = Split("""common/dgraph.c common/dgraphTraversal.c
  common/utils.c common/info.c common/undirPartitioning.c
  common/option.c common/debug.c common/clustering.c
  common/dgraphReader.c common/dgraphDotReader.cpp
  recBisection/vcycle2way.c recBisection/initialBisection.c
  recBisection/rvcycle.c recBisection/dgraphBisection.c recBisection/refinementBis.c
  recBisection/dagP.c""")

dgraphlibsrcs = [os.path.join('#/src', x) for x in dgraphlibsrcs]

l = env.StaticLibrary('lib/'+dgraphlibname, dgraphlibsrcs)

env.Alias('dagplib', l)


rmlgp = env.Program('exe/rMLGP', 'src/recBisection/rMLGP.c')


env.Alias('rMLGP', rmlgp)

# --- Tools  ------------
ginfo = env.Program('exe/ginfo', 'src/common/ginfo.c')

env.Alias('ginfo', ginfo)

env.Alias('tools', [ginfo])


env.Default(rmlgp, ginfo)

env.Alias('all', ['dagplib', 'rMLGP', 'tools'])

env.SConsignFile()
