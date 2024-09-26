cmake_minimum_required (VERSION 3.8)

project(lunar)

add_library(${PROJECT_NAME}
    alt_az.cpp ades2mpc.cpp astfuncs.cpp big_vsop.cpp
    brentmin.cpp cgi_func.cpp classel.cpp conbound.c cospar.cpp date.cpp
    delta_t.cpp de_plan.cpp dist_pa.cpp eart2000.cpp elp82dat.cpp
    eop_prec.cpp getplane.cpp get_time.cpp jsats.cpp lunar2.cpp miscell.cpp moid.cpp
    mpc_code.cpp mpc_fmt.cpp mpc_fmt2.cpp nanosecs.cpp nutation.cpp
    obliquit.cpp pluto.cpp precess.cpp showelem.cpp
    snprintf.cpp sof.cpp spline.cpp ssats.cpp triton.cpp unpack.cpp vislimit.cpp vsopson.cpp
)

add_compile_options(-Wextra -Wall -O3 -pedantic)

target_include_directories(${PROJECT_NAME} PUBLIC . )