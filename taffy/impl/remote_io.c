#include "remote_io.h"
#include "sonLib.h"
#include "htslib/hfile.h"
#include <string.h>
#include <stdlib.h>
#include <errno.h>

bool is_url(const char *path) {
    if (path == NULL) return false;
    return strstr(path, "://") != NULL;
}

char *tai_path_for(const char *input) {
    assert(input != NULL);
    size_t n = strlen(input);
    char *ret = (char *) st_calloc(n + 5, sizeof(char));
    sprintf(ret, "%s.tai", input);
    return ret;
}

/* Slurp an entire URL into a malloc'd buffer. Returns the buffer (caller frees)
 * via *out_buf and sets *out_len. Returns true on success, false on failure
 * (with an error printed to stderr). */
static bool slurp_url(const char *url, char **out_buf, size_t *out_len) {
    hFILE *hf = hopen(url, "r");
    if (hf == NULL) {
        fprintf(stderr, "Unable to open URL: %s (%s)\n", url, strerror(errno));
        fprintf(stderr, "  URL inputs require htslib built with libcurl support.\n"
                        "  If you built htslib yourself, rerun ./configure --enable-libcurl and rebuild.\n");
        return false;
    }
    size_t cap = 64 * 1024;
    size_t len = 0;
    char *buf = (char *) malloc(cap);
    if (buf == NULL) { hclose_abruptly(hf); return false; }
    for (;;) {
        if (len == cap) {
            cap *= 2;
            char *nb = (char *) realloc(buf, cap);
            if (nb == NULL) { free(buf); hclose_abruptly(hf); return false; }
            buf = nb;
        }
        ssize_t n = hread(hf, buf + len, cap - len);
        if (n < 0) {
            fprintf(stderr, "Read error fetching URL: %s\n", url);
            free(buf);
            hclose_abruptly(hf);
            return false;
        }
        if (n == 0) break;
        len += (size_t) n;
    }
    if (hclose(hf) != 0) {
        fprintf(stderr, "Close error after fetching URL: %s\n", url);
        free(buf);
        return false;
    }
    *out_buf = buf;
    *out_len = len;
    return true;
}

FILE *open_tai_for_reading(const char *path) {
    if (!is_url(path)) {
        return fopen(path, "r");
    }
    char *buf = NULL;
    size_t len = 0;
    if (!slurp_url(path, &buf, &len)) {
        return NULL;
    }
    /* Spool to an anonymous tmpfile rather than fmemopen, because downstream
     * LI_construct calls fileno()+bgzf_dopen which needs a real file descriptor.
     * tmpfile() is auto-deleted when closed, so no cleanup required by the caller. */
    FILE *fh = tmpfile();
    if (fh == NULL) {
        fprintf(stderr, "tmpfile failed for buffered .tai from %s\n", path);
        free(buf);
        return NULL;
    }
    if (fwrite(buf, 1, len, fh) != len) {
        fprintf(stderr, "Short write spooling .tai to tmpfile from %s\n", path);
        fclose(fh);
        free(buf);
        return NULL;
    }
    free(buf);
    rewind(fh);
    return fh;
}
