/*
 * sonListPrivate.h
 *
 *  Created on: 28 Jun 2012
 *      Author: benedictpaten
 */

#ifndef SONLISTPRIVATE_H_
#define SONLISTPRIVATE_H_

struct _stList {
    void **list;
    int64_t length;
    int64_t maxLength;
    void (*destructElement)(void *);
};

struct _stListIterator {
    stList *list;
    int64_t index;
};

#endif /* SONLISTPRIVATE_H_ */
