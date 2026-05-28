#ifndef TAF_REMOTE_IO_H_
#define TAF_REMOTE_IO_H_

/*
 * Helpers for reading TAF/MAF + .tai inputs from URLs (HTTP/HTTPS/S3/GCS/etc.)
 * via htslib's hFILE backend. The data file is read on demand through htslib's
 * BGZF/HTTP-Range layer (see line_iterator.c -- LI_construct_from_path goes
 * through bgzf_open which handles both local paths and URLs). The .tai is
 * small and is fetched in full into an in-memory buffer that's then exposed
 * to the rest of the code as a FILE*.
 */

#include <stdio.h>
#include <stdbool.h>

/* Heuristic URL detection: anything containing "://" is treated as a URL.
 * Matches http, https, s3, gs, ftp, etc. -- whatever htslib's hopen supports. */
bool is_url(const char *path);

/* Build the .tai sibling path/URL for a given input. Caller frees the returned
 * string. This is the same construction as tai_path() but with no assumptions
 * about local-vs-remote: it just appends ".tai". */
char *tai_path_for(const char *input);

/* Open an input source for the .tai index, whether local file or URL.
 *
 * For a local path: behaves like fopen(path, "r").
 * For a URL: streams the entire response via htslib's hopen/hread into an
 *            anonymous tmpfile() (rather than fmemopen, since downstream
 *            LI_construct calls fileno()+bgzf_dopen and needs a real fd).
 *            The tmpfile is auto-deleted when the caller fcloses it.
 *
 * Returns NULL on failure (with an error printed to stderr).
 */
FILE *open_tai_for_reading(const char *path);

#endif
