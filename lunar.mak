# Basic astronomical functions library.  Use
# nmake -f lunar.mak [BITS_32=Y] [DLL=Y]
# to build a version that is either 32 or 64 bit;
# and which either builds the library as a DLL or statically

EXES= astcheck.exe astephem.exe calendar.exe cosptest.exe dist.exe \
      easter.exe get_test.exe htc20b.exe jd.exe jevent.exe \
      jpl2b32.exe jsattest.exe lun_test.exe marstime.exe oblitest.exe \
      persian.exe phases.exe prectest.exe ps_1996.exe \
      relativi.exe ssattest.exe tables.exe \
      testprec.exe test_ref.exe uranus1.exe utc_test.exe

all: $(EXES)

LIB_OBJS= ades2mpc.obj alt_az.obj astfuncs.obj \
      big_vsop.obj brentmin.obj classel.obj  \
                   com_file.obj cospar.obj date.obj \
      de_plan.obj delta_t.obj dist_pa.obj  \
      elp82dat.obj eop_prec.obj getplane.obj \
      get_time.obj jsats.obj lunar2.obj  \
      miscell.obj mpc_code.obj mpc_fmt.obj moid.obj \
      nutation.obj obliquit.obj pluto.obj precess.obj  \
      refract.obj refract4.obj rocks.obj showelem.obj sof.obj \
      snprintf.obj \
      spline.obj ssats.obj triton.obj vislimit.obj vsopson.obj

LINK=link /nologo

!ifdef BITS_32
BASE_FLAGS=-nologo -W3 -O2 -MT -D_CRT_SECURE_NO_WARNINGS
LIBNAME=lunar
RM=del
!else
BASE_FLAGS=-nologo -W3 -O2 -MT -D_CRT_SECURE_NO_WARNINGS
LIBNAME=lunar64
RM=del
!endif

!ifdef DLL
COMMON_FLAGS = $(BASE_FLAGS) -EHsc -LD -DNDEBUG
!else
COMMON_FLAGS = $(BASE_FLAGS)
!endif

.cpp.obj:
   cl -c $(COMMON_FLAGS) $<

$(LIBNAME).lib: $(LIB_OBJS)
   del $(LIBNAME).lib
!ifdef DLL
   del $(LIBNAME).dll
   link /DLL /MAP /IMPLIB:$(LIBNAME).lib /DEF:$(LIBNAME).def $(LIB_OBJS)
!else
   lib /OUT:$(LIBNAME).lib $(LIB_OBJS)
!endif

clean:
   $(RM) $(LIB_OBJS)
   $(RM) $(EXES)
   $(RM) ades2mpc.obj astcheck.obj astephem.obj
   $(RM) calendar.obj cosptest.obj dist.obj
   $(RM) easter.obj get_test.obj htc20b.obj jd.obj jevent.obj
   $(RM) jpl2b32.obj jsattest.obj lun_test.obj marstime.obj
   $(RM) mpc_code.obj mpc_fmt.obj oblitest.obj
   $(RM) persian.obj phases.obj ps_1996.obj prectest.obj
   $(RM) relativi.obj ssattest.obj tables.obj
   $(RM) testprec.obj test_ref.obj uranus1.obj utc_test.obj
   $(RM) eart2000.obj gust86.obj lun_tran.obj mpcorb.obj obliqui2.obj
   $(RM) riseset3.obj sof.obj solseqn.obj spline.obj
   $(RM) $(LIBNAME).lib $(LIBNAME).map $(LIBNAME).exp
   $(RM) $(LIBNAME).dll

astcheck.exe: astcheck.obj eart2000.obj mpcorb.obj $(LIBNAME).lib
   $(LINK)    astcheck.obj eart2000.obj mpcorb.obj $(LIBNAME).lib

astephem.exe: astephem.obj eart2000.obj mpcorb.obj $(LIBNAME).lib
   $(LINK)    astephem.obj eart2000.obj mpcorb.obj $(LIBNAME).lib

calendar.exe: calendar.obj $(LIBNAME).lib
   $(LINK)    calendar.obj $(LIBNAME).lib

cosptest.exe: cosptest.obj $(LIBNAME).lib
   $(LINK)    cosptest.obj $(LIBNAME).lib

dist.exe:  dist.obj
   $(LINK) dist.obj

easter.exe: easter.cpp snprintf.obj
   cl -DTEST_CODE $(BASE_FLAGS) easter.cpp snprintf.obj

get_test.exe: get_test.obj $(LIBNAME).lib
   $(LINK)    get_test.obj $(LIBNAME).lib

htc20b.exe: htc20b.cpp
   cl -DTEST_MAIN $(BASE_FLAGS) htc20b.cpp

integrat.exe: integrat.obj $(LIBNAME).lib
   $(LINK)    integrat.obj $(LIBNAME).lib jpleph.lib

jevent.exe: jevent.obj $(LIBNAME).lib
   $(LINK)  jevent.obj $(LIBNAME).lib

jd.exe:    jd.obj $(LIBNAME).lib
   $(LINK) jd.obj $(LIBNAME).lib

jpl2b32.exe: jpl2b32.obj
   $(LINK)   jpl2b32.obj

jsattest.exe: jsattest.obj $(LIBNAME).lib
   $(LINK)    jsattest.obj $(LIBNAME).lib

lun_test.exe: lun_test.obj lun_tran.obj riseset3.obj $(LIBNAME).lib
   $(LINK)    lun_test.obj lun_tran.obj riseset3.obj $(LIBNAME).lib

marstime.exe: marstime.cpp
   cl /DTEST_PROGRAM $(BASE_FLAGS) marstime.cpp

moidtest.exe: moidtest.obj $(LIBNAME).lib
   $(LINK)    moidtest.obj $(LIBNAME).lib

oblitest.exe: oblitest.obj obliqui2.obj spline.obj $(LIBNAME).lib
   $(LINK)    oblitest.obj obliqui2.obj spline.obj $(LIBNAME).lib

persian.exe: persian.obj solseqn.obj $(LIBNAME).lib
   $(LINK)   persian.obj solseqn.obj $(LIBNAME).lib

phases.exe: phases.obj $(LIBNAME).lib
   $(LINK)  phases.obj $(LIBNAME).lib

prectest.exe: prectest.obj $(LIBNAME).lib
   $(LINK)    prectest.obj $(LIBNAME).lib

ps_1996.exe: ps_1996.obj $(LIBNAME).lib
   $(LINK)   ps_1996.obj $(LIBNAME).lib

relativi.exe: relativi.obj $(LIBNAME).lib
   $(LINK)    relativi.obj $(LIBNAME).lib

relativi.obj:
   cl -c $(BASE_FLAGS) /DTEST_CODE relativi.cpp

ssattest.exe: ssattest.obj $(LIBNAME).lib
   $(LINK)    ssattest.obj $(LIBNAME).lib

ssats.obj: ssats.cpp
   cl -c $(COMMON_FLAGS) ssats.cpp

tables.exe: tables.obj riseset3.obj $(LIBNAME).lib
   $(LINK)  tables.obj riseset3.obj $(LIBNAME).lib

testprec.exe: testprec.obj $(LIBNAME).lib
   $(LINK)    testprec.obj $(LIBNAME).lib

test_ref.exe: test_ref.obj refract.obj refract4.obj
   $(LINK)    test_ref.obj refract.obj refract4.obj

uranus1.exe: uranus1.obj gust86.obj
   $(LINK)   uranus1.obj gust86.obj

utc_test.exe: utc_test.obj $(LIBNAME).lib
   $(LINK)    utc_test.obj $(LIBNAME).lib

INSTALL_DIR=..\myincl

install:
   copy afuncs.h   $(INSTALL_DIR)
   copy brentmin.h $(INSTALL_DIR)
   copy cgi_func.h $(INSTALL_DIR)
   copy colors.h   $(INSTALL_DIR)
   copy comets.h   $(INSTALL_DIR)
   copy date.h     $(INSTALL_DIR)
   copy get_bin.h  $(INSTALL_DIR)
   copy gust86.h   $(INSTALL_DIR)
   copy lunar.h    $(INSTALL_DIR)
   copy lun_tran.h $(INSTALL_DIR)
   copy mpc_func.h $(INSTALL_DIR)
   copy riseset3.h $(INSTALL_DIR)
   copy showelem.h $(INSTALL_DIR)
   copy vislimit.h $(INSTALL_DIR)
   copy watdefs.h  $(INSTALL_DIR)
