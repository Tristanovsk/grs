export F2PY=/miniconda3/bin/f2py3.10
# please use an updated version of gcc
# in HAL env:
# module load gcc

DIRS = grs/fortran/grs grs/fortran/grs_a


CLEANDIRS = $(DIRS:%=clean-%)

all: $(DIRS)
$(DIRS):
	$(MAKE) -C $@

clean: $(CLEANDIRS)
$(CLEANDIRS): 
	$(MAKE) -C $(@:clean-%=%) clean

.PHONY: subdirs $(DIRS)
.PHONY: subdirs $(CLEANDIRS)
.PHONY: all clean
