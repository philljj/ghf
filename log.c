#include <stdio.h>
#include "log.h"



void
log_error(const char * errstr)
{
    fprintf(stderr, "error: %s\n", errstr);
    return;
}
