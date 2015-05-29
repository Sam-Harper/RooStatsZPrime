PKG_OBJ_DIR=$(OBJ_DIR)/RooStatsZPrime/src

ROOSTATSZPRIME_LIBFILES        =$(PKG_OBJ_DIR)/DataBox.o $(PKG_OBJ_DIR)/DataPruner.o $(PKG_OBJ_DIR)/ModelConfigurator.o $(PKG_OBJ_DIR)/ModelConfiguratorZprime.o $(PKG_OBJ_DIR)/Pixie.o $(PKG_OBJ_DIR)/PoiRangeEstimator.o $(PKG_OBJ_DIR)/Resultator.o

ROOSTATSZPRIME_LIBNAME 	= $(LIB_DIR)/libRooStatsZPrime.so

$(ROOSTATSZPRIME_LIBNAME):	$(ROOSTATSZPRIME_LIBFILES)
		@$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@
		@echo "$@ done"

RooStatsZPrime/src/%Dict.cc: PACKAGENAME:= RooStatsZPrime
RooStatsZPrime/src/%_LinkDef.h:   packages/RooStatsZPrime/dict/%_LinkDef.h packages/RooStatsZPrime/include/%.hh
#	@echo "joy and happyness $*"
	@sleep 0.001s #dummy command, okay its been a long battle and this hack works... 

STD_LIBS	+= $(ROOSTATSZPRIME_LIBNAME)
