CC 		= gcc
CFLAGS  = -g -pg -Wall
LIBS    = -lm

SRCS = RACIPE.c pcg_basic.c RACIPELIB.c rkf45.c
OBJS = $(SRCS:.c=.o)

MAIN = RACIPE

.PHONY: depend clean

all: $(MAIN)

$(MAIN):$(OBJS)
		$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^
