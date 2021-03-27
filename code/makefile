CC			=				gcc
CFLAGS		=				-Wall -g
LFLAGS		=				-lm

OBJECTS		=				LinkCConfig.o									\
							LinkCTypes.o									\
							LinkCUtilities.o								\
							LinkCHashMap.o								\
							LinkC.o
BINARIES	=				linkc

################

all:	$(BINARIES)

clean:
	/bin/rm -f -r $(BINARIES) $(OBJECTS)

################

LinkCConfig.o:			LinkCConfig.c LinkCConfig.h
	$(CC) $(CFLAGS) $< -c -o $@

LinkCTypes.o:				LinkCTypes.c LinkCTypes.h LinkCConfig.h
	$(CC) $(CFLAGS) $< -c -o $@

LinkCUtilities.o:				LinkCUtilities.c LinkCUtilities.h
	$(CC) $(CFLAGS) $< -c -o $@

LinkCHashMap.o:				LinkCHashMap.c LinkCHashMap.h			\
								LinkCUtilities.h
	$(CC) $(CFLAGS) $< -c -o $@

LinkCMorgan.o:			LinkC.c LinkCConfig.h LinkCTypes.h
	$(CC) $(CFLAGS) $< -c -o $@

linkc:				$(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LFLAGS)
