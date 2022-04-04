# GNU MAKE Makefile for 'lunar' basic astronomical functions library
#  (use 'gmake' for BSD,  and probably 'gmake CLANG=Y')
#
# Usage: make -f [path\]linlunar.mak [CLANG=Y] [W64=Y] [W32=Y] [MSWIN=Y] [tgt]
#
# where tgt can be any of:
# [all|astcheck|astephem|calendar... clean]
# [install|install_integrat|install_astcheck]
#
#	'W64'/'W32' = cross-compile for Windows 64 or 32 bits,  using MinGW,
#    on a Linux or BSD box
#	'MSWIN' = compile for Windows,  using MinGW,  on a Windows machine
#	'CLANG' = use clang instead of GCC;  BSD/Linux only
# None of these: compile using g++ on BSD or Linux
#	Note that I've only tried clang on PC-BSD (which is based on FreeBSD).
#
# 'integrat' is not built as part of 'make'.  If you want that,  run
# 'make integrat' and then,  optionally,  'make install_integrat'.

# As CC is an implicit variable, a simple CC?=gcc doesn't work.
# We have to use this trick from https://stackoverflow.com/a/42958970
ifeq ($(origin CPP),default)
	CPP=g++
endif
ifeq ($(origin CC),default)
	CC=gcc
endif
LIBSADDED=
EXE=
CFLAGS+=-Wextra -Wall -O3 -pedantic

ifdef DEBUG
	CFLAGS += -g
endif

ifdef ERRORS
	CFLAGS += -Werror
endif

# You can have your include files in ~/include and libraries in
# ~/lib,  in which case only the current user can use them;  or
# (with root privileges) you can install them to /usr/local/include
# and /usr/local/lib for all to enjoy.

PREFIX?=~
ifdef GLOBAL
	INSTALL_DIR=/usr/local
else
	INSTALL_DIR=$(PREFIX)
endif

ifdef CLANG
	CPP=clang++
	CC=clang
endif

RM=rm -f

ifeq ($(OS),Windows_NT)
    detected_OS := Windows
else
    detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

ifeq ($(detected_OS),Linux)
    CP = cp -u
else
    CP = cp
endif


# I'm using 'mkdir -p' to avoid error messages if the directory exists.
# It may fail on very old systems,  and will probably fail on non-POSIX
# systems.  If so,  change to '-mkdir' and ignore errors.

ifdef MSWIN
	EXE=.exe
	MKDIR=-mkdir
else
	LIBSADDED=-lm
	MKDIR=mkdir -p
endif

LIB_DIR=$(INSTALL_DIR)/lib

ifdef W64
   CPP=x86_64-w64-mingw32-g++
   CC=x86_64-w64-mingw32-gcc
   EXE=.exe
   LIB_DIR=$(INSTALL_DIR)/win_lib
   LIBSADDED=-L $(LIB_DIR) -mwindows
endif

ifdef W32
   CPP=i686-w64-mingw32-g++
   CC=i686-w64-mingw32-gcc
   EXE=.exe
   LIB_DIR=$(INSTALL_DIR)/win_lib32
   LIBSADDED=-L $(LIB_DIR) -mwindows
endif

ifeq ($(SHARED),Y)
	LIBEXE = $(CC)
	CFLAGS += -fPIC
	LIBFLAGS = -shared -Wl,-soname,liblunar.so.1.0.1 -lc -o
	LIBLUNAR = liblunar.so.1.0.1
else
	LIBEXE = ar
	LIBFLAGS = crsv
	LIBLUNAR = liblunar.a
endif

all: add_off$(EXE) adestest$(EXE) astcheck$(EXE) astephem$(EXE) \
   calendar$(EXE) cgicheck$(EXE) chinese$(EXE) colors$(EXE) \
   colors2$(EXE) cosptest$(EXE) csv2ades$(EXE) dist$(EXE) \
   easter$(EXE) get_test$(EXE) gtest$(EXE) htc20b$(EXE) jd$(EXE)\
   jevent$(EXE) jpl2b32$(EXE) jsattest$(EXE) lun_test$(EXE) \
   marstime$(EXE) moidtest$(EXE) mpc2sof$(EXE) oblitest$(EXE) \
   persian$(EXE) parallax$(EXE) parallax.cgi phases$(EXE) \
   prectest$(EXE) prectes2$(EXE) ps_1996$(EXE) ssattest$(EXE) \
   tables$(EXE) test_des$(EXE) test_ref$(EXE) testprec$(EXE) \
   themis$(EXE) them_cat$(EXE) uranus1$(EXE) utc_test$(EXE)

install:
	$(MKDIR) $(INSTALL_DIR)/include
	$(CP) afuncs.h   $(INSTALL_DIR)/include
	$(CP) brentmin.h $(INSTALL_DIR)/include
	$(CP) cgi_func.h $(INSTALL_DIR)/include
	$(CP) comets.h   $(INSTALL_DIR)/include
	$(CP) date.h     $(INSTALL_DIR)/include
	$(CP) get_bin.h  $(INSTALL_DIR)/include
	$(CP) lunar.h    $(INSTALL_DIR)/include
	$(CP) mpc_func.h $(INSTALL_DIR)/include
	$(CP) showelem.h $(INSTALL_DIR)/include
	$(CP) vislimit.h $(INSTALL_DIR)/include
	$(CP) watdefs.h  $(INSTALL_DIR)/include
	$(MKDIR) $(LIB_DIR)
	$(CP) $(LIBLUNAR) $(LIB_DIR)
	$(MKDIR) $(INSTALL_DIR)/bin

install_astcheck:
	$(CP) astcheck$(EXE) $(INSTALL_DIR)/bin

install_integrat:
	$(CP) integrat $(INSTALL_DIR)/bin

uninstall:
	rm -f $(INSTALL_DIR)/include/afuncs.h
	rm -f $(INSTALL_DIR)/include/cgi_func.h
	rm -f $(INSTALL_DIR)/include/comets.h
	rm -f $(INSTALL_DIR)/include/date.h
	rm -f $(INSTALL_DIR)/include/get_bin.h
	rm -f $(INSTALL_DIR)/include/lunar.h
	rm -f $(INSTALL_DIR)/include/mpc_func.h
	rm -f $(INSTALL_DIR)/include/showelem.h
	rm -f $(INSTALL_DIR)/include/vislimit.h
	rm -f $(INSTALL_DIR)/include/watdefs.h
	rm -f $(INSTALL_DIR)/lib/$(LIBLUNAR)
	rm -f $(INSTALL_DIR)/bin/astcheck$(EXE)
	rm -f $(INSTALL_DIR)/bin/integrat$(EXE)

.cpp.o:
	$(CPP) $(CFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

OBJS= alt_az.o ades2mpc.o astfuncs.o big_vsop.o  \
   brentmin.o cgi_func.o classel.o cospar.o date.o  \
   delta_t.o de_plan.o dist_pa.o eart2000.o elp82dat.o \
   eop_prec.o getplane.o get_time.o jsats.o lunar2.o miscell.o moid.o \
   mpc_code.o mpc_fmt.o nutation.o obliquit.o pluto.o precess.o showelem.o \
   sof.o spline.o ssats.o triton.o unpack.o vislimit.o vsopson.o

$(LIBLUNAR): $(OBJS)
	$(LIBEXE) $(LIBFLAGS) $(LIBLUNAR) $(OBJS)

clean:
	$(RM) $(OBJS)
	$(RM) adestest.o add_off.o astcheck.o astephem.o calendar.o cgicheck.o
	$(RM) cosptest.o csv2ades.o get_test.o gtest.o gust86.o htc20b.o integrat.o jd.o
	$(RM) jevent.o jpl2b32.o jsattest.o lun_test.o lun_tran.o mms.o
	$(RM) moidtest.o mpcorb.o oblitest.o obliqui2.o persian.o phases.o
	$(RM) prectes2.o prectest.o ps_1996.o refract.o refract4.o riseset3.o solseqn.o
	$(RM) ssattest.o tables.o test_des.o test_ref.o testprec.o
	$(RM) themis.o transit.o uranus1.o utc_test.o
	$(RM) add_off$(EXE) add_off.cgi
	$(RM) adestest$(EXE) astcheck$(EXE) astephem$(EXE) calendar$(EXE)
	$(RM) cgicheck$(EXE) chinese$(EXE) colors$(EXE)
	$(RM) colors2$(EXE) cosptest$(EXE) csv2ades$(EXE) dist$(EXE)
	$(RM) easter$(EXE) get_test$(EXE) gtest$(EXE) htc20b$(EXE)
	$(RM) integrat$(EXE) jd$(EXE) jevent$(EXE) jpl2b32$(EXE)
	$(RM) jsattest$(EXE) lun_test$(EXE) marstime$(EXE) moidtest$(EXE) mms$(EXE)
	$(RM) mpc2sof$(EXE) oblitest$(EXE) parallax$(EXE) parallax.cgi
	$(RM) persian$(EXE) phases$(EXE) prectest$(EXE) prectes2$(EXE)
	$(RM) ps_1996$(EXE) relativi$(EXE) solseqn$(EXE) ssattest$(EXE) tables$(EXE)
	$(RM) test_des$(EXE) test_ref$(EXE) testprec$(EXE) themis$(EXE)
	$(RM) them_cat$(EXE) transit$(EXE) uranus1$(EXE) utc_test$(EXE) $(LIBLUNAR)

add_off$(EXE): add_off.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o add_off$(EXE) add_off.o $(LIBLUNAR) $(LIBSADDED)

add_off.cgi: add_off.c $(LIBLUNAR)
	$(CC) $(CFLAGS) -o add_off.cgi -DON_LINE_VERSION add_off.c $(LIBLUNAR) $(LIBSADDED)

adestest$(EXE): adestest.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o adestest$(EXE) adestest.o $(LIBLUNAR) $(LIBSADDED)

astcheck$(EXE): astcheck.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o astcheck$(EXE) astcheck.o $(LIBLUNAR) $(LIBSADDED)

astephem$(EXE): astephem.o mpcorb.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o astephem$(EXE) astephem.o mpcorb.o $(LIBLUNAR) $(LIBSADDED)

calendar$(EXE): calendar.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o calendar$(EXE) calendar.o   $(LIBLUNAR) $(LIBSADDED)

cgicheck$(EXE): astcheck.cpp $(LIBLUNAR) cgicheck.o
	$(CC) $(CFLAGS) -o cgicheck$(EXE) -DCGI_VERSION cgicheck.o astcheck.cpp $(LIBLUNAR) $(LIBSADDED)

chinese$(EXE): chinese.cpp
	$(CC) $(CFLAGS) -o chinese$(EXE) chinese.cpp

colors$(EXE): colors.cpp
	$(CC) $(CFLAGS) -o colors$(EXE) colors.cpp -DSIMPLE_TEST_PROGRAM

colors2$(EXE): colors2.cpp
	$(CC) $(CFLAGS) -o colors2$(EXE) colors2.cpp -DTEST_FUNC

cosptest$(EXE): cosptest.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o cosptest$(EXE) cosptest.o   $(LIBLUNAR) $(LIBSADDED)

csv2ades$(EXE): csv2ades.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o csv2ades$(EXE) csv2ades.o   $(LIBLUNAR) $(LIBSADDED)

dist$(EXE): dist.cpp
	$(CC) $(CFLAGS) -o dist$(EXE) dist.cpp $(LIBSADDED)

easter$(EXE): easter.cpp $(LIBLUNAR)
	$(CC) $(CFLAGS) -o easter$(EXE) -DTEST_CODE easter.cpp $(LIBLUNAR) $(LIBSADDED)

get_test$(EXE): get_test.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o get_test$(EXE) get_test.o $(LIBLUNAR) $(LIBSADDED)

gtest$(EXE): gtest.c
	$(CC) $(FLAGS) -o gtest$(EXE) gtest.c $(LIBSADDED)

htc20b$(EXE): htc20b.cpp $(LIBLUNAR)
	$(CC) $(CFLAGS) -o htc20b$(EXE) -DTEST_MAIN htc20b.cpp $(LIBLUNAR) $(LIBSADDED)

integrat$(EXE): integrat.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o integrat$(EXE) integrat.o $(LIBLUNAR) $(LIBSADDED) -L $(INSTALL_DIR)/lib -ljpl

integrat.o: integrat.cpp
	$(CC) $(CFLAGS) -c -I $(INSTALL_DIR)/include $<

jd$(EXE): jd.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o jd$(EXE) jd.o $(LIBLUNAR) $(LIBSADDED)

jevent$(EXE):                    jevent.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o jevent$(EXE) jevent.o $(LIBLUNAR) $(LIBSADDED)

jpl2b32$(EXE):                    jpl2b32.o
	$(CC) $(CFLAGS) -o jpl2b32$(EXE) jpl2b32.o $(LIBSADDED)

jsattest$(EXE): jsattest.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o jsattest$(EXE) jsattest.o $(LIBLUNAR) $(LIBSADDED)

lun_test$(EXE): lun_test.o lun_tran.o riseset3.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o lun_test$(EXE) lun_test.o lun_tran.o riseset3.o $(LIBLUNAR) $(LIBSADDED)

marstime$(EXE): marstime.cpp
	$(CC) $(CFLAGS) -o marstime$(EXE) marstime.cpp -DTEST_PROGRAM $(LIBSADDED)

mms$(EXE):                    mms.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o mms$(EXE) mms.o $(LIBLUNAR) $(LIBSADDED)

moidtest$(EXE): moidtest.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o moidtest$(EXE) moidtest.o $(LIBLUNAR) $(LIBSADDED)

mpc2sof$(EXE): mpc2sof.cpp mpcorb.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o mpc2sof$(EXE) mpc2sof.cpp mpcorb.o $(LIBLUNAR) $(LIBSADDED)

mpc_code$(EXE): mpc_code.cpp
	$(CC) $(CFLAGS) -o mpc_code$(EXE) mpc_code.cpp -DTEST_CODE

oblitest$(EXE): oblitest.o obliqui2.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o oblitest$(EXE) oblitest.o obliqui2.o $(LIBLUNAR) $(LIBSADDED)

parallax.cgi: parallax.cpp $(LIBLUNAR)
	$(CC) $(CFLAGS) -o parallax.cgi parallax.cpp $(LIBLUNAR) $(LIBSADDED) -DCGI_VERSION

parallax$(EXE): parallax.cpp $(LIBLUNAR)
	$(CC) $(CFLAGS) -o parallax$(EXE) parallax.cpp $(LIBLUNAR) $(LIBSADDED)

persian$(EXE): persian.o solseqn.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o persian$(EXE) persian.o solseqn.o $(LIBLUNAR) $(LIBSADDED)

phases$(EXE): phases.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o phases$(EXE)   phases.o   $(LIBLUNAR) $(LIBSADDED)

prectest$(EXE): prectest.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o prectest$(EXE) prectest.o $(LIBLUNAR) $(LIBSADDED)

prectes2$(EXE): prectes2.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o prectes2$(EXE) prectes2.o $(LIBLUNAR) $(LIBSADDED)

ps_1996$(EXE): ps_1996.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o ps_1996$(EXE)   ps_1996.o   $(LIBLUNAR) $(LIBSADDED)

relativi$(EXE): relativi.cpp $(LIBLUNAR)
	$(CC) $(CFLAGS) -o relativi$(EXE) -DTEST_CODE relativi.cpp $(LIBLUNAR) $(LIBSADDED)

sof$(EXE): sof.cpp $(LIBLUNAR)
	$(CC) $(CFLAGS) -DTEST_CODE -o sof$(EXE) sof.cpp -lm $(LIBLUNAR)

spline$(EXE): spline.cpp
	$(CC) $(CFLAGS) -DTEST_CODE -o spline$(EXE) spline.cpp -lm

ssattest$(EXE): ssattest.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o ssattest$(EXE) ssattest.o $(LIBLUNAR) $(LIBSADDED)

tables$(EXE):                    tables.o riseset3.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o tables$(EXE) tables.o riseset3.o $(LIBLUNAR) $(LIBSADDED)

test_des$(EXE):                    test_des.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o test_des$(EXE) test_des.o $(LIBLUNAR) $(LIBSADDED)

test_ref$(EXE):                    test_ref.o refract.o refract4.o
	$(CC) $(CFLAGS) -o test_ref$(EXE) test_ref.o refract.o refract4.o $(LIBSADDED)

testprec$(EXE):                    testprec.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o testprec$(EXE) testprec.o $(LIBLUNAR) $(LIBSADDED)

test_min$(EXE):                    test_min.o brentmin.o
	$(CC) $(CFLAGS) -o test_min$(EXE) test_min.o brentmin.o $(LIBSADDED)

them_cat$(EXE): them_cat.c
	$(CC) $(FLAGS) -o them_cat$(EXE) them_cat.c

themis$(EXE):                    themis.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o themis$(EXE) themis.o $(LIBLUNAR) $(LIBSADDED)

transit$(EXE):                    transit.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o transit$(EXE) transit.o $(LIBLUNAR) $(LIBSADDED)

uranus1$(EXE): uranus1.o gust86.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o uranus1$(EXE) uranus1.o gust86.o $(LIBSADDED)

utc_test$(EXE):                utc_test.o $(LIBLUNAR)
	$(CC) $(CFLAGS) -o utc_test$(EXE) utc_test.o $(LIBLUNAR) $(LIBSADDED)

