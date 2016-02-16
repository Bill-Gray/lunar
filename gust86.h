/* gust86.h: header file for Uranian satellite computations
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

/* Output state vector is in AU and AU/second */

#ifndef __stdcall
#define __stdcall
#endif

#ifdef __cplusplus
extern "C" {
#endif
void __stdcall gust86_posn( const double jde, const int isat, double *r );
#ifdef __cplusplus
}
#endif

#define GUST86_ARIEL          0
#define GUST86_UMBRIEL        1
#define GUST86_TITANIA        2
#define GUST86_OBERON         3
#define GUST86_MIRANDA        4
