# GNU MAKE Makefile for 'lunar' basic astronomical functions library
#  (use 'gmake' for BSD,  and probably 'gmake CLANG=Y')
#
# Usage: make -f [path\]linlunar.mak [CLANG=Y] [XCOMPILE=Y] [MSWIN=Y] [tgt]
#
# where tgt can be any of:
# [all|astcheck|astephem|calendar... clean]
#
#	'XCOMPILE' = cross-compile for Windows,  using MinGW,  on a Linux or BSD box
#	'MSWIN' = compile for Windows,  using MinGW,  on a Windows machine
#	'CLANG' = use clang instead of GCC;  BSD/Linux only
# None of these: compile using g++ on BSD or Linux
#	Note that I've only tried clang on PC-BSD (which is based on FreeBSD).

CC=g++
LIBSADDED=
EXE=
CFLAGS=-Wextra -Wall -O3 -pedantic -Wno-unused-parameter

ifdef CLANG
	CC=clang
endif

RM=rm -f

ifdef MSWIN
	EXE=.exe
else
	LIBSADDED=-lm
endif

ifdef XCOMPILE
   CC=x86_64-w64-mingw32-g++
   EXE=.exe
   LIBSADDED=
endif

all: astcheck$(EXE) astephem$(EXE) calendar$(EXE) cgicheck$(EXE)  \
   colors$(EXE) colors2$(EXE) cosptest$(EXE) dist$(EXE) easter$(EXE) \
   get_test$(EXE) htc20b$(EXE) jd$(EXE) \
   jevent$(EXE) jpl2b32$(EXE) jsattest$(EXE) lun_test$(EXE) \
   marstime$(EXE) oblitest$(EXE) persian$(EXE) phases$(EXE) \
   prectest$(EXE) ps_1996$(EXE) ssattest$(EXE) tables$(EXE) \
   test_ref$(EXE) testprec$(EXE) uranus1$(EXE) utc_test$(EXE)

install:
	cp afuncs.h   /usr/local/include
	cp comets.h   /usr/local/include
	cp showelem.h /usr/local/include
	cp date.h     /usr/local/include
	cp lunar.h    /usr/local/include
	cp watdefs.h  /usr/local/include
	cp liblunar.a /usr/local/lib

uninstall:
	rm -f /usr/local/include/afuncs.h
	rm -f /usr/local/include/comets.h
	rm -f /usr/local/include/showelem.h
	rm -f /usr/local/include/date.h
	rm -f /usr/local/include/watdefs.h
	rm -f /usr/local/lib/liblunar.a

.cpp.o:
	$(CC) $(CFLAGS) -c $<

OBJS= alt_az.o astfuncs.o big_vsop.o classel.o cospar.o date.o delta_t.o \
   de_plan.o dist_pa.o eart2000.o elp82dat.o \
   eop_prec.o getplane.o get_time.o \
   jsats.o lunar2.o miscell.o nutation.o obliquit.o pluto.o precess.o \
   showelem.o spline.o ssats.o triton.o vislimit.o vsopson.o

liblunar.a: $(OBJS)
	ar crsv liblunar.a $(OBJS)

clean:
	$(RM) $(OBJS)
	$(RM) astcheck.o astephem.o calendar.o cgicheck.o cgi_func.o
	$(RM) cosptest.o get_test.o gust86.o htc20b.o integrat.o jd.o
	$(RM) jevent.o jpl2b32.o jsattest.o lun_test.o lun_tran.o
	$(RM) mpcorb.o oblitest.o obliqui2.o persian.o phases.o
	$(RM) prectest.o ps_1996.o refract.o refract4.o riseset3.o solseqn.o
	$(RM) ssattest.o tables.o test_ref.o testprec.o uranus1.o utc_test.o
	$(RM) astcheck$(EXE) astephem$(EXE) calendar$(EXE) cgicheck$(EXE) colors$(EXE)
	$(RM) colors2$(EXE) cosptest$(EXE) dist$(EXE) easter$(EXE) get_test$(EXE)
	$(RM) htc20b$(EXE) integrat$(EXE) jd$(EXE) jevent$(EXE) jpl2b32$(EXE) jsattest$(EXE)
	$(RM) lun_test$(EXE) marstime$(EXE) oblitest$(EXE) persian$(EXE) phases$(EXE) prectest$(EXE)
	$(RM) ps_1996$(EXE) relativi$(EXE) solseqn$(EXE) ssattest$(EXE) tables$(EXE)
	$(RM) test_ref$(EXE) testprec$(EXE) uranus1$(EXE) utc_test$(EXE) liblunar.a

astcheck$(EXE): astcheck.o mpcorb.o liblunar.a
	$(CC) $(CFLAGS) -o astcheck$(EXE) astcheck.o mpcorb.o liblunar.a $(LIBSADDED)

astephem$(EXE): astephem.o mpcorb.o liblunar.a
	$(CC) $(CFLAGS) -o astephem$(EXE) astephem.o mpcorb.o liblunar.a $(LIBSADDED)

calendar$(EXE): calendar.o liblunar.a
	$(CC) $(CFLAGS) -o calendar$(EXE) calendar.o   liblunar.a $(LIBSADDED)

cgicheck$(EXE): astcheck.cpp mpcorb.o liblunar.a cgicheck.o cgi_func.o
	$(CC) $(CFLAGS) -o cgicheck$(EXE) -DCGI_VERSION cgicheck.o astcheck.cpp mpcorb.o cgi_func.o liblunar.a $(LIBSADDED)

colors$(EXE): colors.cpp
	$(CC) $(CFLAGS) -o colors$(EXE) colors.cpp -DSIMPLE_TEST_PROGRAM

colors2$(EXE): colors2.cpp
	$(CC) $(CFLAGS) -o colors2$(EXE) colors2.cpp -DTEST_FUNC

cosptest$(EXE): cosptest.o liblunar.a
	$(CC) $(CFLAGS) -o cosptest$(EXE) cosptest.o   liblunar.a $(LIBSADDED)

dist$(EXE): dist.cpp
	$(CC) $(CFLAGS) -o dist$(EXE) dist.cpp $(LIBSADDED)

easter$(EXE): easter.cpp liblunar.a
	$(CC) $(CFLAGS) -o easter$(EXE) -DTEST_CODE easter.cpp liblunar.a $(LIBSADDED)

get_test$(EXE): get_test.o liblunar.a
	$(CC) $(CFLAGS) -o get_test$(EXE) get_test.o liblunar.a $(LIBSADDED)

htc20b$(EXE): htc20b.cpp liblunar.a
	$(CC) $(CFLAGS) -o htc20b$(EXE) -DTEST_MAIN htc20b.cpp liblunar.a $(LIBSADDED)

integrat$(EXE): integrat.o liblunar.a
	$(CC) $(CFLAGS) -o integrat$(EXE) integrat.o liblunar.a $(LIBSADDED) -ljpl

jd$(EXE): jd.o liblunar.a
	$(CC) $(CFLAGS) -o jd$(EXE) jd.o liblunar.a $(LIBSADDED)

jevent$(EXE):                    jevent.o liblunar.a
	$(CC) $(CFLAGS) -o jevent$(EXE) jevent.o liblunar.a $(LIBSADDED)

jpl2b32$(EXE):                    jpl2b32.o
	$(CC) $(CFLAGS) -o jpl2b32$(EXE) jpl2b32.o $(LIBSADDED)

jsattest$(EXE): jsattest.o liblunar.a
	$(CC) $(CFLAGS) -o jsattest$(EXE) jsattest.o liblunar.a $(LIBSADDED)

lun_test$(EXE): lun_test.o lun_tran.o riseset3.o liblunar.a
	$(CC) $(CFLAGS) -o lun_test$(EXE) lun_test.o lun_tran.o riseset3.o liblunar.a $(LIBSADDED)

marstime$(EXE): marstime.cpp
	$(CC) $(CFLAGS) -o marstime$(EXE) marstime.cpp -DTEST_PROGRAM $(LIBSADDED)

oblitest$(EXE): oblitest.o obliqui2.o liblunar.a
	$(CC) $(CFLAGS) -o oblitest$(EXE) oblitest.o obliqui2.o liblunar.a $(LIBSADDED)

persian$(EXE): persian.o solseqn.o liblunar.a
	$(CC) $(CFLAGS) -o persian$(EXE) persian.o solseqn.o liblunar.a $(LIBSADDED)

phases$(EXE): phases.o liblunar.a
	$(CC) $(CFLAGS) -o phases$(EXE)   phases.o   liblunar.a $(LIBSADDED)

prectest$(EXE): prectest.o liblunar.a
	$(CC) $(CFLAGS) -o prectest$(EXE) prectest.o liblunar.a $(LIBSADDED)

ps_1996$(EXE): ps_1996.o liblunar.a
	$(CC) $(CFLAGS) -o ps_1996$(EXE)   ps_1996.o   liblunar.a $(LIBSADDED)

relativi$(EXE): relativi.cpp liblunar.a
	$(CC) $(CFLAGS) -o relativi$(EXE) -DTEST_CODE relativi.cpp liblunar.a $(LIBSADDED)

ssattest$(EXE): ssattest.o liblunar.a
	$(CC) $(CFLAGS) -o ssattest$(EXE) ssattest.o liblunar.a $(LIBSADDED)

tables$(EXE):                    tables.o riseset3.o liblunar.a
	$(CC) $(CFLAGS) -o tables$(EXE) tables.o riseset3.o liblunar.a $(LIBSADDED)

test_ref$(EXE):                    test_ref.o refract.o refract4.o
	$(CC) $(CFLAGS) -o test_ref$(EXE) test_ref.o refract.o refract4.o $(LIBSADDED)

testprec$(EXE):                    testprec.o liblunar.a
	$(CC) $(CFLAGS) -o testprec$(EXE) testprec.o liblunar.a $(LIBSADDED)

uranus1$(EXE): uranus1.o gust86.o
	$(CC) $(CFLAGS) -o uranus1$(EXE) uranus1.o gust86.o $(LIBSADDED)

utc_test$(EXE):                utc_test.o liblunar.a
	$(CC) $(CFLAGS) -o utc_test$(EXE) utc_test.o liblunar.a $(LIBSADDED)

