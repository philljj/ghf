#if !defined(OPTIONS_H)
#define OPTIONS_H

#include <stdbool.h>
#include <stdlib.h>

#include "spin_types.h"

enum hf_type_t {
    RHF,
    UHF,
    GHF
};

typedef enum hf_type_t hf_type_t;

struct options_t {
    bool      debug;
    bool      quiet;
    bool      in_memory;
    size_t    threads;
    size_t    sleep_time;
    hf_type_t type;
    char *    input_file;
};

bool      validate_options(int argc, char * argv[]);
bool      is_debug(void);
bool      is_quiet(void);
bool      in_memory(void);
size_t    num_threads(void);
size_t    get_sleep_time(void);
hf_type_t get_hf_type(void);
char *    get_input_file(void);

#endif
