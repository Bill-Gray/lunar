/* mini_dll.cpp: probably obsolete code for a 16-bit DLL!
I don't think use has been made of this in a decade or so.

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

// libentry.cpp : Defines the entry point for the DLL application.
//
/* Usually,  the following lines come from a separate file,  STDAFX.H.
   However,  for the 'mini-DLL',  STDAFX is used only in this file,
   and maintaining a separate include file seems superfluous.       */

/*        ------ STDAFX.H BEGINS HERE ------ */
// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__298970FE_2CDF_11D5_A298_000102CA0D3E__INCLUDED_)
#define AFX_STDAFX_H__298970FE_2CDF_11D5_A298_000102CA0D3E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


// Insert your headers here
#define WIN32_LEAN_AND_MEAN      // Exclude rarely-used stuff from Windows headers

#include <windows.h>

// TODO: reference additional headers your program requires here

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__298970FE_2CDF_11D5_A298_000102CA0D3E__INCLUDED_)
/*        ------ STDAFX.H ENDS HERE ------ */

BOOL APIENTRY DllMain( HANDLE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                )
{
    return TRUE;
}

