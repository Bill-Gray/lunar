# Basic astronomical functions library - Win32 .DLL version

EXES= astcheck.exe astephem.exe calendar.exe cosptest.exe dist.exe \
      easter.exe get_test.exe htc20b.exe jd.exe jevent.exe \
      jpl2b32.exe jsattest.exe lun_test.exe marstime.exe oblitest.exe \
      persian.exe phases.exe ps_1996.exe relativi.exe ssattest.exe tables.exe \
      testprec.exe test_ref.exe uranus1.exe utc_test.exe

all: $(EXES)

LIB_OBJS= alt_az.obj astfuncs.obj big_vsop.obj classel.obj com_file.obj \
      cospar.obj date.obj de_plan.obj delta_t.obj dist_pa.obj elp82dat.obj \
      getplane.obj get_time.obj jsats.obj lunar2.obj  \
      miscell.obj nutation.obj obliquit.obj pluto.obj precess.obj  \
      refract.obj refract4.obj rocks.obj showelem.obj \
      ssats.obj triton.obj vislimit.obj vsopson.obj

LINK=link /nologo

!ifdef BITS_32
COMMON_FLAGS=-nologo -W3
LIBNAME=lunar
RM=rm
!else
COMMON_FLAGS=-nologo -W3 -D_CRT_SECURE_NO_WARNINGS
LIBNAME=lunar64
RM=del
!endif

$(LIBNAME).lib: $(LIB_OBJS)
   del $(LIBNAME).lib
   del $(LIBNAME).dll
   link /DLL /MAP /IMPLIB:$(LIBNAME).lib /DEF:$(LIBNAME).def $(LIB_OBJS)

clean:
   $(RM) $(LIB_OBJS)
   $(RM) $(EXES)
   $(RM) astcheck.obj astephem.obj calendar.obj cosptest.obj dist.obj
   $(RM) easter.obj get_test.obj htc20b.obj jd.obj jevent.obj
   $(RM) jpl2b32.obj jsattest.obj lun_test.obj marstime.obj oblitest.obj
   $(RM) persian.obj phases.obj ps_1996.obj relativi.obj ssattest.obj tables.obj
   $(RM) testprec.obj test_ref.obj uranus1.obj utc_test.obj
   $(RM) eart2000.obj gust86.obj lun_tran.obj mpcorb.obj obliqui2.obj
   $(RM) riseset3.obj solseqn.obj spline.obj
   $(RM) $(LIBNAME).lib $(LIBNAME).map $(LIBNAME).exp

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

easter.exe: easter.cpp
   cl -Ox -DTEST_CODE $(COMMON_FLAGS) easter.cpp

get_test.exe: get_test.obj $(LIBNAME).lib
   $(LINK)    get_test.obj $(LIBNAME).lib

htc20b.exe: htc20b.cpp
   cl -Ox -DTEST_MAIN $(COMMON_FLAGS) -nologo htc20b.cpp

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
   cl /Ox /DTEST_PROGRAM $(COMMON_FLAGS) marstime.cpp

oblitest.exe: oblitest.obj obliqui2.obj spline.obj $(LIBNAME).lib
   $(LINK)    oblitest.obj obliqui2.obj spline.obj $(LIBNAME).lib

persian.exe: persian.obj solseqn.obj $(LIBNAME).lib
   $(LINK)   persian.obj solseqn.obj $(LIBNAME).lib

phases.exe: phases.obj $(LIBNAME).lib
   $(LINK)  phases.obj $(LIBNAME).lib

ps_1996.exe: ps_1996.obj $(LIBNAME).lib
   $(LINK)   ps_1996.obj $(LIBNAME).lib

relativi.exe: relativi.obj $(LIBNAME).lib
   $(LINK)    relativi.obj $(LIBNAME).lib

relativi.obj:
   cl /c /Od /DTEST_CODE $(COMMON_FLAGS) relativi.cpp

ssattest.exe: ssattest.obj $(LIBNAME).lib
   $(LINK)    ssattest.obj $(LIBNAME).lib

ssats.obj: ssats.cpp
   cl -Od -EHsc -c -LD $(COMMON_FLAGS) ssats.cpp

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

# CFLAGS=-Ox -GX -c -LD  $(COMMON_FLAGS)
CFLAGS=-Ox -EHsc -c -LD -DNDEBUG $(COMMON_FLAGS)

.cpp.obj:
   cl $(CFLAGS) $<

