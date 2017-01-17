SHELL = /bin/csh
CC = gcc
#CFLAGS = -I. # general operation
CFLAGS = -I. -g # basic debug
#CFLAGS = -I. -g -Wall # fuller debug
LIBRARY = -lm -lz
INSTALL_DIR = /depot/phig/apps/linux

all: 
	make VicOutputASMStats
	make GetVicHeader
	make VicBinaryDump2ASCII
	make create_xmask
	make compute_average_grid_forcings

clean::
	/bin/rm -f *.o core *.core log *~

VicOutputASMStats: VicOutputASMStats.c DateFuncs.c VicUtilities.c VicStats_1d.c
	$(CC) -o $(INSTALL_DIR)/VicOutputASMStats VicOutputASMStats.c DateFuncs.c VicUtilities.c VicStats_1d.c $(CFLAGS) $(LIBRARY)

GetVicHeader: GetVicHeader.c VicUtilities.c 
	$(CC) -o $(INSTALL_DIR)/GetVicHeader GetVicHeader.c VicUtilities.c $(CFLAGS) $(LIBRARY)

VicBinaryDump2ASCII: VicBinaryDump2ASCII.c VicUtilities.c
	$(CC) -o $(INSTALL_DIR)/VicBinaryDump2ASCII VicBinaryDump2ASCII.c VicUtilities.c $(CFLAGS) $(LIBRARY)

create_xmask: create_xmask.c
	$(CC) -o $(INSTALL_DIR)/create_xmask create_xmask.c $(CFLAGS) $(LIBRARY)

compute_average_grid_forcings: compute_average_grid_forcings.c
	$(CC) -o $(INSTALL_DIR)/compute_average_grid_forcings compute_average_grid_forcings.c $(CFLAGS) $(LIBRARY)

