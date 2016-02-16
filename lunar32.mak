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

lunar.lib: $(LIB_OBJS)
   del lunar.lib
   lib /OUT:lunar.lib $(LIB_OBJS)

clean:
   rm $(LIB_OBJS)
   rm $(EXES)
   rm astcheck.obj astephem.obj calendar.obj cosptest.obj dist.obj
   rm easter.obj get_test.obj htc20b.obj jd.obj jevent.obj
   rm jpl2b32.obj jsattest.obj lun_test.obj marstime.obj oblitest.obj
   rm persian.obj phases.obj ps_1996.obj relativi.obj ssattest.obj tables.obj
   rm testprec.obj test_ref.obj uranus1.obj utc_test.obj
   rm eart2000.obj gust86.obj lun_tran.obj mpcorb.obj obliqui2.obj
   rm riseset3.obj solseqn.obj spline.obj
   rm lunar.lib lunar.map

astcheck.exe: astcheck.obj eart2000.obj mpcorb.obj lunar.lib
   $(LINK)    astcheck.obj eart2000.obj mpcorb.obj lunar.lib

astephem.exe: astephem.obj eart2000.obj mpcorb.obj lunar.lib
   $(LINK)    astephem.obj eart2000.obj mpcorb.obj lunar.lib

calendar.exe: calendar.obj lunar.lib
   $(LINK)    calendar.obj lunar.lib

cosptest.exe: cosptest.obj lunar.lib
   $(LINK)    cosptest.obj lunar.lib

dist.exe:  dist.obj
   $(LINK) dist.obj

easter.exe: easter.cpp
   cl -W4 -Ox -DTEST_CODE -nologo easter.cpp

get_test.exe: get_test.obj lunar.lib
   $(LINK)    get_test.obj lunar.lib

htc20b.exe: htc20b.cpp
   cl -W3 -Ox -DTEST_MAIN -nologo htc20b.cpp

jevent.exe: jevent.obj lunar.lib
   $(LINK)  jevent.obj lunar.lib

jd.exe:    jd.obj lunar.lib
   $(LINK) jd.obj lunar.lib

jpl2b32.exe: jpl2b32.obj
   $(LINK)   jpl2b32.obj

jsattest.exe: jsattest.obj lunar.lib
   $(LINK)    jsattest.obj lunar.lib

lun_test.exe: lun_test.obj lun_tran.obj riseset3.obj lunar.lib
   $(LINK)    lun_test.obj lun_tran.obj riseset3.obj lunar.lib

marstime.exe: marstime.cpp
   cl /nologo /Ox /W3 /DTEST_PROGRAM marstime.cpp

oblitest.exe: oblitest.obj obliqui2.obj spline.obj lunar.lib
   $(LINK)    oblitest.obj obliqui2.obj spline.obj lunar.lib

persian.exe: persian.obj solseqn.obj lunar.lib
   $(LINK)   persian.obj solseqn.obj lunar.lib

phases.exe: phases.obj lunar.lib
   $(LINK)  phases.obj lunar.lib

ps_1996.exe: ps_1996.obj lunar.lib
   $(LINK)   ps_1996.obj lunar.lib

relativi.exe: relativi.obj lunar.lib
   $(LINK)    relativi.obj lunar.lib

relativi.obj:
   cl /c /Od /W3 /DTEST_CODE -nologo relativi.cpp

ssattest.exe: ssattest.obj lunar.lib
   $(LINK)    ssattest.obj lunar.lib

ssats.obj: ssats.cpp
   cl -W3 -Od -GX -c -LD -nologo -I\myincl ssats.cpp

tables.exe: tables.obj riseset3.obj lunar.lib
   $(LINK)  tables.obj riseset3.obj lunar.lib

testprec.exe: testprec.obj lunar.lib
   $(LINK)    testprec.obj lunar.lib

test_ref.exe: test_ref.obj refract.obj refract4.obj
   $(LINK)    test_ref.obj refract.obj refract4.obj

uranus1.exe: uranus1.obj gust86.obj
   $(LINK)   uranus1.obj gust86.obj

utc_test.exe: utc_test.obj lunar.lib
   $(LINK)    utc_test.obj lunar.lib

CFLAGS=-W3 -Ox -GX -c -LD  -nologo -D_CRT_SECURE_NO_WARNINGS

.cpp.obj:
   cl $(CFLAGS) $<

