/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.3 2003/07/25 13:52:00 karypis Exp $
 */

/*
#define	DEBUG		1
#define	DMALLOC		1
*/

/* added by fkummer ****************************************************************/

// surpress "deprecated"-warning
#pragma warning (disable: 4996) 
// surpress "unreference local variable"-warning
#pragma warning (disable: 4101)
// surpress "conversion from ... to ..., possible loss of data"-warning
#pragma warning (disable: 4244)
// surpress "conversion from ... to ..., possible loss of data"-warning
#pragma warning (disable: 4305)

/***************************************************************** end fkummer ****/


#include "stdheaders.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include "../ParMETIS/parmetis.h"  /* Get the idxtype definition */
#include "defs.h"
#include "struct.h"
#include "macros.h"
#include "rename.h"
#include "proto.h"

