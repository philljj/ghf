#include <getopt.h>
#include <stdio.h>

#include "options.h"

#define GHF_MAX_THREADS 64

static void set_default_options(void);

static struct options_t opts;

bool
validate_options(int    argc,
                 char * argv[])
{
    set_default_options();

    int opt = 0;

    while ((opt = getopt(argc, argv, "dmqsf:t:?")) != -1) {
        switch (opt) {
        case 'd':
            opts.debug = true;
            break;

        case 'm':
            opts.in_memory = true;
            break;

        case 'f':
            opts.input_file = optarg;
            break;

        case 't':
            opts.threads = strtoul(optarg, 0, 10);
            printf("using threads: %zu\n", opts.threads);
            break;

        case 'q':
            opts.quiet = true;
            break;

        case 's':
            opts.sleep_time = 1;
            printf("using sleep: %zu\n", opts.sleep_time);
            break;

        case '?':
        default:
            return false;
        }
    }

    if (num_threads() == 0) {
        fprintf(stderr, "error: threads == 0 is invalid\n");
        return false;
    }

    if (num_threads() > GHF_MAX_THREADS ) {
        fprintf(stderr, "error: threads > %d is invalid\n", GHF_MAX_THREADS);
        return false;
    }

    if (!opts.quiet) {
        printf("using input file: %s\n", get_input_file());
    }

    return true;
}



static void
set_default_options(void)
{
    opts.debug = false;
    opts.quiet = false;
    opts.in_memory = false;
    opts.threads = 1;
    opts.sleep_time = 0;
    opts.type = RHF;
    opts.input_file = 0;

    return;
}



bool
is_debug(void)
{
    return opts.debug;
}



bool
is_quiet(void)
{
    return opts.quiet;
}



bool
in_memory(void)
{
    return opts.in_memory;
}



size_t
num_threads(void)
{
    return opts.threads;
}



size_t
get_sleep_time(void)
{
    return opts.sleep_time;
}



enum hf_type
get_hf_type(void)
{
    return opts.type;
}



char *
get_input_file(void)
{
    return opts.input_file;
}
