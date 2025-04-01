CC = gcc
CFLAGS = -Wall -O2 -std=c99 -frounding-math -ffp-contract=off
LDFLAGS = -lm


TARGET = test_convert

SRCS = test_convert.c convert.c
OBJS = $(SRCS:.c=.o)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean