#include <sys/queue.h>

// AI-aided implementation of a bottom-up iterative mergesort macro for
// TAILQ lists, inspired by the works of Simon Tatham
#define TAILQ_MERGESORT(head, headname, type, field, cmp) do {         \
    struct type *p, *q, *e, *tail;                                     \
    int insize = 1, nmerges, psize, qsize, i;                          \
    if (TAILQ_EMPTY(head) || TAILQ_NEXT(TAILQ_FIRST(head), field) == NULL) \
        break;                                                         \
                                                                       \
    while (1) {                                                        \
        p = TAILQ_FIRST(head);                                         \
        TAILQ_INIT(head);                                              \
        tail = NULL;                                                   \
        nmerges = 0;                                                   \
                                                                       \
        while (p) {                                                    \
            nmerges++;                                                 \
            q = p;                                                     \
            psize = 0;                                                 \
            for (i = 0; i < insize; i++) {                             \
                psize++;                                               \
                q = TAILQ_NEXT(q, field);                              \
                if (!q) break;                                         \
            }                                                          \
            qsize = insize;                                            \
                                                                       \
            while (psize > 0 || (qsize > 0 && q)) {                    \
                if (psize == 0) {                                      \
                    e = q; q = TAILQ_NEXT(q, field); qsize--;          \
                } else if (qsize == 0 || !q) {                         \
                    e = p; p = TAILQ_NEXT(p, field); psize--;          \
                } else if (cmp(p, q) <= 0) {                           \
                    e = p; p = TAILQ_NEXT(p, field); psize--;          \
                } else {                                               \
                    e = q; q = TAILQ_NEXT(q, field); qsize--;          \
                }                                                      \
                                                                       \
                if (tail) {                                            \
                    TAILQ_INSERT_AFTER(head, tail, e, field);          \
                } else {                                               \
                    TAILQ_INSERT_HEAD(head, e, field);                 \
                }                                                      \
                tail = e;                                              \
            }                                                          \
            p = q;                                                     \
        }                                                              \
        if (nmerges <= 1) break;                                       \
        insize *= 2;                                                   \
    }                                                                  \
} while(0)
