# Basic astronomical functions library.  Use
# nmake -f lunar.mak [BITS_32=Y] [DLL=Y]
# to build a version that is either 32 or 64 bit;
# and which either builds the library as a DLL or statically

EXES= add_off.exe adestest.exe astcheck.exe astephem.exe \
      calendar.exe chinese.exe colors.exe colors2.exe cosptest.exe csv2ades.exe dist.exe \
      easter.exe get_test.exe gtest.exe htc20b.exe jd.exe jevent.exe \
      jpl2b32.exe jsattest.exe lun_test.exe marstime.exe \
      moidtest.exe mpc_time.exe mpc2sof.exe oblitest.exe parallax.exe \
      persian.exe phases.exe prectest.exe prectes2.exe ps_1996.exe \
      relativi.exe ssattest.exe tables.exe test_des.exe test_ref.exe \
      testprec.exe test_ref.exe themis.exe them_cat.exe uranus1.exe utc_test.exe

all: $(EXES)

LIB_OBJS= ades2mpc.obj alt_az.obj astfuncs.obj \
      big_vsop.obj brentmin.obj classel.obj  \
      com_file.obj conbound.obj cospar.obj date.obj \
      de_plan.obj delta_t.obj dist_pa.obj  \
      elp82dat.obj eop_prec.obj getplane.obj \
      get_time.obj jsats.obj lunar2.obj miscell.obj mpc_code.obj \
      mpc_fmt.obj mpc_fmt2.obj moid.obj nanosecs.obj \
      nutation.obj obliquit.obj pluto.obj precess.obj  \
      refract.obj refract4.obj rocks.obj showelem.obj sof.obj \
      snprintf.obj spline.obj ssats.obj \
      unpack.obj triton.obj vislimit.obj vsopson.obj

LINK=link /nologo

!ifdef BITS_32
BASE_FLAGS=-nologo -W3 -O2 -MT -D_CRT_SECURE_NO_WARNINGS
LIBNAME=lunar
JPLLIBNAME=jpleph32
RM=del
!else
BASE_FLAGS=-nologo -W3 -O2 -MT -D_CRT_SECURE_NO_WARNINGS
LIBNAME=lunar64
JPLLIBNAME=jpleph64
RM=del
!endif

!ifdef DLL
COMMON_FLAGS = $(BASE_FLAGS) -EHsc -LD -DNDEBUG
!else
COMMON_FLAGS = $(BASE_FLAGS)
!endif

.cpp.obj:
   cl -c $(COMMON_FLAGS) $<

.c.obj:
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
   $(RM) add_off.obj ades2mpc.obj adestest.obj astcheck.obj
   $(RM) astephem.obj calendar.obj chinese.obj colors.obj
   $(RM) colors2.obj csv2ades.obj cosptest.obj dist.obj
   $(RM) eart2000.obj easter.obj get_test.obj gtest.obj
   $(RM) gust86.obj htc20b.obj jd.obj jevent.obj
   $(RM) jpl2b32.obj jsattest.obj lun_test.obj lun_tran.obj
   $(RM) marstime.obj mpc_code.obj mpc_fmt.obj mpc_fmt2.obj mpcorb.obj
   $(RM) moidtest.obj mpc_time.obj mpc2sof.obj obliqui2.obj oblitest.obj
   $(RM) parallax.obj persian.obj phases.obj ps_1996.obj
   $(RM) prectest.obj prectes2.obj relativi.obj riseset3.obj
   $(RM) sof.obj solseqn.obj spline.obj ssattest.obj tables.obj
   $(RM) testprec.obj test_des.obj test_ref.obj themis.obj
   $(RM) them_cat.obj uranus1.obj utc_test.obj
   $(RM) $(LIBNAME).lib $(LIBNAME).map $(LIBNAME).exp
   $(RM) $(LIBNAME).dll

add_off.exe:  add_off.obj $(LIBNAME).lib
   $(LINK)    add_off.obj $(LIBNAME).lib urlmon.lib

adestest.exe: adestest.obj $(LIBNAME).lib
   $(LINK)    adestest.obj $(LIBNAME).lib

astcheck.exe: astcheck.obj eart2000.obj mpcorb.obj $(LIBNAME).lib
   $(LINK)    astcheck.obj eart2000.obj mpcorb.obj $(LIBNAME).lib

astephem.exe: astephem.obj eart2000.obj mpcorb.obj $(LIBNAME).lib
   $(LINK)    astephem.obj eart2000.obj mpcorb.obj $(LIBNAME).lib

calendar.exe: calendar.obj $(LIBNAME).lib
   $(LINK)    calendar.obj $(LIBNAME).lib

chinese.exe: chinese.cpp snprintf.obj
   cl -DTEST_CODE $(BASE_FLAGS) chinese.cpp snprintf.obj

colors.exe: colors.cpp
   cl -DSIMPLE_TEST_PROGRAM $(BASE_FLAGS) colors.cpp

colors2.exe: colors2.cpp
   cl -DTEST_FUNC $(BASE_FLAGS) colors2.cpp

cosptest.exe: cosptest.obj $(LIBNAME).lib
   $(LINK)    cosptest.obj $(LIBNAME).lib

csv2ades.exe: csv2ades.obj $(LIBNAME).lib
   $(LINK)    csv2ades.obj $(LIBNAME).lib

dist.exe:  dist.obj
   $(LINK) dist.obj

easter.exe: easter.cpp snprintf.obj
   cl -DTEST_CODE $(BASE_FLAGS) easter.cpp snprintf.obj

get_test.exe: get_test.obj $(LIBNAME).lib
   $(LINK)    get_test.obj $(LIBNAME).lib

htc20b.exe: htc20b.cpp
   cl -DTEST_MAIN $(BASE_FLAGS) htc20b.cpp

integrat.exe: integrat.obj $(LIBNAME).lib
   $(LINK)    integrat.obj $(LIBNAME).lib $(JPLLIBNAME).lib

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

marstime.exe: marstime.cpp snprintf.obj
   cl /DTEST_PROGRAM $(BASE_FLAGS) marstime.cpp snprintf.obj

moidtest.exe: moidtest.obj $(LIBNAME).lib
   $(LINK)    moidtest.obj $(LIBNAME).lib

mpc_time.exe: mpc_time.obj $(LIBNAME).lib
   $(LINK)    mpc_time.obj $(LIBNAME).lib

mpc2sof.exe: mpc2sof.obj mpcorb.obj $(LIBNAME).lib
   $(LINK)   mpc2sof.obj mpcorb.obj $(LIBNAME).lib

oblitest.exe: oblitest.obj obliqui2.obj spline.obj $(LIBNAME).lib
   $(LINK)    oblitest.obj obliqui2.obj spline.obj $(LIBNAME).lib

parallax.exe: parallax.obj $(LIBNAME).lib
   $(LINK)    parallax.obj $(LIBNAME).lib

persian.exe: persian.obj solseqn.obj $(LIBNAME).lib
   $(LINK)   persian.obj solseqn.obj $(LIBNAME).lib

phases.exe: phases.obj $(LIBNAME).lib
   $(LINK)  phases.obj $(LIBNAME).lib

prectest.exe: prectest.obj $(LIBNAME).lib
   $(LINK)    prectest.obj $(LIBNAME).lib

prectes2.exe: prectes2.obj $(LIBNAME).lib
   $(LINK)    prectes2.obj $(LIBNAME).lib

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

test_des.exe: test_des.obj $(LIBNAME).lib
   $(LINK)    test_des.obj $(LIBNAME).lib

test_ref.exe: test_ref.obj refract.obj refract4.obj
   $(LINK)    test_ref.obj refract.obj refract4.obj

them_cat.exe: them_cat.c snprintf.obj
   cl -DTEST_CODE $(BASE_FLAGS) them_cat.c snprintf.obj

themis.exe:   themis.obj $(LIBNAME).lib
   $(LINK)    themis.obj $(LIBNAME).lib

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
   copy stringex.h $(INSTALL_DIR)
   copy vislimit.h $(INSTALL_DIR)
   copy watdefs.h  $(INSTALL_DIR)
