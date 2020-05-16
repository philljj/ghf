#if !defined(OPTIONS_H)
#define OPTIONS_H

#include <stdbool.h>
#include <stdlib.h>

enum hf_type {
    RHF,
    UHF,
    GHF
};

struct options_t {
    bool         debug;
    bool         quiet;
    bool         in_memory;
    size_t       threads;
    size_t       sleep_time;
    enum hf_type type;
    char *       input_file;
};

bool         validate_options(int argc, char * argv[]);
bool         is_debug(void);
bool         is_quiet(void);
bool         in_memory(void);
size_t       num_threads(void);
size_t       get_sleep_time(void);
enum hf_type get_hf_type(void);
char *       get_input_file(void);

#endif
