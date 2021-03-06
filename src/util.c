#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"



void *
safer_calloc(size_t       count,
             size_t       size,
             const char * what)
{
    if (!count || !size) {
        fprintf(stderr, "fatal: invalid calloc(%zu, %zu)\n", count, size);
        exit(1);
    }

    void * p = calloc(count, size);

    if (!p) {
        if (what && *what) {
            fprintf(stderr, "fatal: %s: calloc(%zu, %zu) failed\n", what,
                    count, size);
        }
        else {
            fprintf(stderr, "fatal: calloc(%zu, %zu) failed\n", count, size);
        }

        exit(1);
    }

    return p;
}
