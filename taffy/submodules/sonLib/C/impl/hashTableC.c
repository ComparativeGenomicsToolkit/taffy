
/* Copyright (C) 2004 Christopher Clark <firstname.lastname@cl.cam.ac.uk> */

#include "hashTableC.h"
#include "hashTablePrivateC.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "sonLibGlobalsInternal.h"

/*
 Credit for primes table: Aaron Krowne
 http://br.endernet.org/~akrowne/
 http://planetmath.org/encyclopedia/GoodHashTablePrimes.html
 */
static const uint64_t primes[] = { 53, 97, 193, 389, 769, 1543, 3079, 6151,
        12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739,
        6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189,
        805306457, 1610612741 };
const uint64_t prime_table_length = sizeof(primes) / sizeof(primes[0]);
const float max_load_factor = 0.65;

/* indexFor */
static uint64_t indexFor(uint64_t tablelength, uint64_t hashvalue) {
    return (hashvalue % tablelength);
}

/*****************************************************************************/
struct hashtable *
create_hashtable(uint64_t minsize, uint64_t(*hashf)(const void*), int(*eqf)(
        const void*, const void*), void(*keyFree)(void *), void(*valueFree)(
        void *)) {
    struct hashtable *h;
    uint64_t pindex, size = primes[0];
    /* Check requested hashtable isn't too large */
    if (minsize > (1u << 30))
        return NULL;
    /* Enforce size as prime */
    for (pindex = 0; pindex < prime_table_length; pindex++) {
        if (primes[pindex] > minsize) {
            size = primes[pindex];
            break;
        }
    }
    h = (struct hashtable *) st_malloc(sizeof(struct hashtable));
    if (NULL == h)
        return NULL; /*oom*/
    h->table = (struct entry **) st_malloc(sizeof(struct entry*) * size);
    if (NULL == h->table) {
        free(h);
        return NULL;
    } /*oom*/
    memset(h->table, 0, size * sizeof(struct entry *));
    h->tablelength = size;
    h->primeindex = pindex;
    h->entrycount = 0;
    h->hashfn = hashf;
    h->eqfn = eqf;
    h->loadlimit = (uint64_t) ceil(size * max_load_factor);
    //my additions
    h->keyFree = keyFree;
    h->valueFree = valueFree;
    return h;
}

/*****************************************************************************/
uint64_t hashP(struct hashtable *h, void *k) {
    /* Aim to protect against poor hash functions by adding logic here
     * - logic taken from java 1.4 hashtable source */
    uint64_t i = h->hashfn(k);
    i += ~(i << 9);
    i ^= ((i >> 14) | (i << 18)); /* >>> */
    i += (i << 4);
    i ^= ((i >> 10) | (i << 22)); /* >>> */
    return i;
}

/*****************************************************************************/
static int64_t hashtable_expand(struct hashtable *h) {
    /* Double the size of the table to accomodate more entries */
    struct entry **newtable;
    struct entry *e;
    //struct entry **pE;
    uint64_t newsize, i, index;
    /* Check we're not hitting max capacity */
    if (h->primeindex == (prime_table_length - 1))
        return 0;
    newsize = primes[++(h->primeindex)];
    newtable = (struct entry **) st_malloc(sizeof(struct entry*) * newsize);
    if (NULL != newtable) {
        memset(newtable, 0, newsize * sizeof(struct entry *));
        /* This algorithm is not 'stable'. ie. it reverses the list
         * when it transfers entries between the tables */
        for (i = 0; i < h->tablelength; i++) {
            while (NULL != (e = h->table[i])) {
                h->table[i] = e->next;
                index = indexFor(newsize, e->h);
                e->next = newtable[index];
                newtable[index] = e;
            }
        }
        free(h->table);
        h->table = newtable;
    }
    /* Plan B: realloc instead */
    else {
        return 0; //the realloc doesn't always work
        /*newtable = (struct entry **)
         realloc(h->table, newsize * sizeof(struct entry *));
         if (NULL == newtable) { (h->primeindex)--; return 0; }
         h->table = newtable;
         memset(newtable[h->tablelength], 0, newsize - h->tablelength);
         for (i = 0; i < h->tablelength; i++) {
         for (pE = &(newtable[i]), e = *pE; e != NULL; e = *pE) {
         index = indexFor(newsize,e->h);
         if (index == i)
         {
         pE = &(e->next);
         }
         else
         {
         *pE = e->next;
         e->next = newtable[index];
         newtable[index] = e;
         }
         }
         } */
    }
    h->tablelength = newsize;
    h->loadlimit = (uint64_t) ceil(newsize * max_load_factor);
    return -1;
}

/*****************************************************************************/
uint64_t hashtable_count(struct hashtable *h) {
    return h->entrycount;
}

/*****************************************************************************/
int64_t hashtable_insert(struct hashtable *h, void *k, void *v) {
    /* This method allows duplicate keys - but they shouldn't be used */
    uint64_t index;
    struct entry *e;
    if (++(h->entrycount) > h->loadlimit) {
        /* Ignore the return value. If expand fails, we should
         * still try cramming just this value INT_32o the existing table
         * -- we may not have memory for a larger table, but one more
         * element may be ok. Next time we insert, we'll try expanding again.*/
        if (!hashtable_expand(h)) {
            fprintf(stderr, "Expand hash failed \n");
            exit(1);
        }
    }
    e = (struct entry *) st_malloc(sizeof(struct entry));
    if (NULL == e) {
        --(h->entrycount);
        return 0;
    } /*oom*/
    e->h = hashP(h, k);
    index = indexFor(h->tablelength, e->h);
    e->k = k;
    e->v = v;
    e->next = h->table[index];
    h->table[index] = e;
    return -1;
}

/*****************************************************************************/
void * /* returns value associated with key */
hashtable_search(struct hashtable *h, void *k) {
    struct entry *e;
    uint64_t hashvalue, index;
    hashvalue = hashP(h, k);
    index = indexFor(h->tablelength, hashvalue);
    e = h->table[index];
    while (NULL != e) {
        /* Check hash value to short circuit heavier comparison */
        if ((hashvalue == e->h) && (h->eqfn(k, e->k)))
            return e->v;
        e = e->next;
    }
    return NULL;
}

/*****************************************************************************/
void * /* returns value associated with key */
hashtable_remove(struct hashtable *h, void *k, int64_t freeKey) {
    /* TODO: consider compacting the table when the load factor drops enough,
     *       or provide a 'compact' method. */

    struct entry *e;
    struct entry **pE;
    void *v;
    uint64_t hashvalue, index;

    hashvalue = hashP(h, k);
    index = indexFor(h->tablelength, hashP(h, k));
    pE = &(h->table[index]);
    e = *pE;
    while (NULL != e) {
        /* Check hash value to short circuit heavier comparison */
        if ((hashvalue == e->h) && (h->eqfn(k, e->k))) {
            *pE = e->next;
            h->entrycount--;
            v = e->v;
            if (freeKey) {
                h->keyFree(e->k);
            }
            free(e);
            return v;
        }
        pE = &(e->next);
        e = e->next;
    }
    return NULL;
}

/*****************************************************************************/
/* destroy */
void hashtable_destroy(struct hashtable *h, int64_t free_values,
        int64_t free_keys) {
    uint64_t i;
    struct entry *e, *f;
    struct entry **table = h->table;
    if (free_keys) {
        if (free_values) {
            for (i = 0; i < h->tablelength; i++) {
                e = table[i];
                while (NULL != e) {
                    f = e;
                    e = e->next;
                    h->keyFree(f->k);
                    h->valueFree(f->v);
                    free(f);
                }
            }
        } else {
            for (i = 0; i < h->tablelength; i++) {
                e = table[i];
                while (NULL != e) {
                    f = e;
                    e = e->next;
                    h->keyFree(f->k);
                    free(f);
                }
            }
        }
    } else {
        if (free_values) {
            for (i = 0; i < h->tablelength; i++) {
                e = table[i];
                while (NULL != e) {
                    f = e;
                    e = e->next;
                    h->valueFree(f->v);
                    free(f);
                }
            }
        } else {
            for (i = 0; i < h->tablelength; i++) {
                e = table[i];
                while (NULL != e) {
                    f = e;
                    e = e->next;
                    free(f);
                }
            }
        }
    }
    free(h->table);
    free(h);
}

/*
 * Copyright (c) 2002, Christopher Clark
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the original author; nor the names of any contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
