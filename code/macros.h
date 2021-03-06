/*====================================================================
    macros.h
 
 	macros used to access data structures and perform quick tests.

  ====================================================================*/

#include <stdlib.h>

/* general-purpose macros */
#define SWAP(t,x,y)	{ t = x; x = y; y = t; }

#define SQR(x) ((x)*(x))

#define NEW(p,type)	if ((p=(type *) malloc (sizeof(type))) == NULL) {\
				printf ("Out of Memory!\n");\
				exit(0);\
			}

#define FREE(p)		if (p) { free ((char *) p); p = NULL; }


//add p in front of head
#define ADD( head, p )  if ( head )  { \
				p->next = head; \
				p->prev = head->prev; \
				head->prev = p; \
				p->prev->next = p; \
			} \
			else { \
				head = p; \
				head->next = head->prev = p; \
			}

//
//delete p
//
//if there is only one element, head, delete head
//if p is head, update head
//
#define DELETE( head, p ) if ( head )  { \
				if ( head == head->next ) \
					head = NULL;  \
				else if ( p == head ) \
					head = head->next; \
				p->next->prev = p->prev;  \
				p->prev->next = p->next;  \
				FREE( p ); \
			} 

