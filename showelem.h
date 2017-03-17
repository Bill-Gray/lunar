/* showelem.h: header file for 8-line element display functions
Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

#define SHOWELEM_PRECISION_MASK      0x0f
#define SHOWELEM_PERIH_TIME_MASK     0x10
#define SHOWELEM_OMIT_PQ_MASK        0x20
#define SHOWELEM_COMET_MAGS_NUCLEAR  0x40

/* REMEMBER:  set 'central_obj', 'epoch', 'abs_mag', 'slope_param' fields */

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

int DLL_FUNC elements_in_mpc_format( char *obuff, const ELEMENTS *elem,
                const char *obj_id, const int is_cometary, const int format);

#ifdef __cplusplus
}
#endif  /* #ifdef __cplusplus */
