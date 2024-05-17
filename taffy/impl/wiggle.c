#include "taf.h"
#include "sonLib.h"

/*
 * Add a tag to the set of tags
 */
static void add_tag(stHash *tags, char *key, char *value) {
    assert(stHash_search(tags, key) == NULL);  // Check we don't have multiple values per key
    stHash_insert(tags, stString_copy(key), stString_copy(value));
}

/*
 * Parses a wiggle header line returning a hash of keys to values, all keys/values are strings.
 */
static stHash *parse_header(char *line) {
    // Hash to hold key:value pairs
    stHash *tags = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);

    // Tokenize the line
    stList *tokens = stString_split(line);

    // Parse the step
    char *step = stList_get(tokens, 0);
    if(strcmp(step, "fixedStep") != 0 && strcmp(step, "variableStep") != 0) {
        st_errAbort("Misformed wiggle header line: %s\n", line);
    }
    add_tag(tags, "fixed_step", strcmp(step, "fixedStep") == 0 ? "1" : "0");

    // Parse the remaining tokens
    for(int64_t i=1; i<stList_length(tokens); i++) {
        char *tag = stList_get(tokens, i);
        stList *tokenized_tag = stString_splitByString(tag, "=");
        if(stList_length(tokenized_tag) != 2) {
            st_errAbort("Misformed wiggle header line tag: %s\n", tag);
        }
        add_tag(tags, stList_get(tokenized_tag, 0), stList_get(tokenized_tag, 1));  // Add
        stList_destruct(tokenized_tag);  // Cleanup
    }

    // Cleanup
    stList_destruct(tokens);

    return tags;
}

/*
 * Get a tag from the set of tags
 */
static char *get_tag(stHash *tags, char *key, char *default_value) {
    char *c = stHash_search(tags, key);
    return c == NULL ? default_value : c;
}

/*
 * Add a coordinate-value pair to values
 */
static void add_coordinate_value(stHash *values, int64_t coordinate, double value) {
    assert(stHash_search(values, (void *)coordinate) == NULL);  // Check we don't have multiple values per coordinate
    double *j = st_malloc(sizeof(double));
    *j = value;
    stHash_insert(values, (void *)coordinate, (void *)j);
}

double wig_get_value(stHash *wig, char *seq, int64_t coordinate, double default_value) {
    stHash *values = stHash_search(wig, seq);
    if(values == NULL) {
        return default_value;
    }
    double *value = stHash_search(values, (void *)coordinate);
    return value == NULL ? default_value : *value;
}

stHash *wig_parse(char *file, char *seq_prefix, bool make_zero_based) {
    // Make a hash of string names to nested hashes, each of which is a hash of int64_t to floats
    stHash *seq_intervals = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                              free, (void (*)(void *))stHash_destruct);

    // Get file handle
    FILE *f = fopen(file, "r");
    if(f == NULL) {
        st_errAbort("Failed to open wig file: %s\n", file);
    }
    LI *li = LI_construct(f);

    char *line = LI_get_next_line(li);
    while (line) {  // While we haven't got to the end of the file
        // First line must be a header
        stHash *header = parse_header(line);
        free(line);

        // Get the sequence name, prepending the seq_prefix
        char *seq_name = stString_print("%s%s", seq_prefix, stHash_search(header, "chrom"));

        // Get the values associated with the sequence name
        stHash *values;
        if ((values = stHash_search(seq_intervals, seq_name)) == NULL) {
            values = stHash_construct2(NULL, free);
            stHash_insert(seq_intervals, stString_copy(seq_name), values);
        }

        // Get the span and type of step (fixed or variable)
        int64_t span = atol(get_tag(header, "span", "1"));
        bool fixed_step = atol(get_tag(header, "fixed_step", "1"));

        if (fixed_step) {  // If a fixed step
            int64_t step = atol(get_tag(header, "step", "1"));
            assert(span <= step);  // The span must be less than or equal to the step or we'd have overlapping values
            assert(stHash_search(header, "start") != NULL);  // For a fixed step there must be a start value
            int64_t i = atol(stHash_search(header, "start")) + (make_zero_based ? -1 : 0);
            assert(i >= 0); // Sanity check for start coordinate
            while (1) {
                line = LI_get_next_line(li);

                // If we've hit the end of the file, stop
                if(line == NULL) {
                    break;
                }

                stList *tokens = stString_split(line); // Tokenize
                if (stList_length(tokens) == 0) { // Case we have a stray empty line
                    stList_destruct(tokens);
                    free(line);
                    continue; // Go back and get the next line
                }
                if (stList_length(tokens) > 1) { // Case we have hit another header line
                    stList_destruct(tokens);
                    break; // Jump to outer loop
                }

                // Case we have an entry
                double f = atof(stList_get(tokens, 0)); // The value
                for (int64_t j = 0; j < span; j++) {
                    add_coordinate_value(values, i + j, f);
                }
                i += step;

                // Clean up
                free(line);
                stList_destruct(tokens);
            }
        } else {
            while (1) {
                line = LI_get_next_line(li);

                // If we've hit the end of the file, stop
                if(line == NULL) {
                    break;
                }

                stList *tokens = stString_split(line);
                if (stList_length(tokens) == 0) { // Case we have a stray empty line
                    free(line);
                    stList_destruct(tokens);
                    continue;
                }
                if (stList_length(tokens) > 2 ||
                    strcmp(stList_get(tokens, 0), "variableStep") == 0 ||
                    strcmp(stList_get(tokens, 0), "fixedStep") == 0) { // Case we have hit another header line
                    stList_destruct(tokens);
                    break;
                }

                // Case we have a value
                int64_t i = atol(stList_get(tokens, 0)) + (make_zero_based ? -1 : 0);
                assert(i >= 0); // Sanity check for coordinate
                double f = atof(stList_get(tokens, 1));
                for (int64_t j = 0; j < span; j++) {
                    add_coordinate_value(values, i + j, f);
                }

                // Cleanup
                free(line);
                stList_destruct(tokens);
            }
        }

        // Clean up
        free(seq_name);
        stHash_destruct(header);
    }
    return seq_intervals;
}

