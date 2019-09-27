## Use parallel BLAS for all runs
export OMP_NUM_THREADS=8
## dd:hh:mm:ss
export MAX_DURATION=05:00:00:00
export MAX_MEMORY=20gb

## Phony
all: publish

## Auto generate file with all targets
generated.mk: Makefile
	rm -rf generated.mk
	echo 'cat("commercial-A1: ", paste0("results_WBScod_Age/cod-commercial-m", 1:8, "_A1", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("commercial-A2: ", paste0("results_WBScod_Age/cod-commercial-m", 1:8, "_A2", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("commercial-A3: ", paste0("results_WBScod_Age/cod-commercial-m", 1:8, "_A3", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("commercial-A4: ", paste0("results_WBScod_Age/cod-commercial-m", 1:8, "_A4", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("commercial-A5: ", paste0("results_WBScod_Age/cod-commercial-m", 1:8, "_A5", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("commercial-A6: ", paste0("results_WBScod_Age/cod-commercial-m", 1:8, "_A6", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("survey-A1: ", paste0("results_WBScod_Age/cod-survey-m", 1:8, "_A1", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("survey-A2: ", paste0("results_WBScod_Age/cod-survey-m", 1:8, "_A2", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("survey-A3: ", paste0("results_WBScod_Age/cod-survey-m", 1:8, "_A3", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("survey-A4: ", paste0("results_WBScod_Age/cod-survey-m", 1:8, "_A4", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("survey-A5: ", paste0("results_WBScod_Age/cod-survey-m", 1:8, "_A5", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("survey-A6: ", paste0("results_WBScod_Age/cod-survey-m", 1:8, "_A6", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("both-A1: ", paste0("results_WBScod_Age/cod-both-m", 1:8, "_A1", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("both-A2: ", paste0("results_WBScod_Age/cod-both-m", 1:8, "_A2", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("both-A3: ", paste0("results_WBScod_Age/cod-both-m", 1:8, "_A3", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("both-A4: ", paste0("results_WBScod_Age/cod-both-m", 1:8, "_A4", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("both-A5: ", paste0("results_WBScod_Age/cod-both-m", 1:8, "_A5", ".RData"))' | R --slave >> $@
	echo >> $@
	echo 'cat("both-A6: ", paste0("results_WBScod_Age/cod-both-m", 1:8, "_A6", ".RData"))' | R --slave >> $@
	echo >> $@

-include generated.mk


## NOTE: if more years have to be included, than include those years for each case in the lines below!
commercial: commercial-A1 commercial-A2 commercial-A3 commercial-A4 commercial-A5 commercial-A6
survey: survey-A1 survey-A2 survey-A3 survey-A4 survey-A5 survey-A6
both: both-A1 both-A2 both-A3 both-A4 both-A5 both-A6

results_WBScod_Age: commercial survey both

## Run the model with inputs
RUN_MODEL = R --slave < LGNB_Rmodel.R

## Build the model
MODEL_SRC = LGNB.cpp
MODEL_BIN = LGNB.so
$(MODEL_BIN) : $(MODEL_SRC)
	echo "TMB:::compile('$(MODEL_SRC)', '-Ofast')" | R --slave

## Target 1
results_WBScod_Age/cod-commercial-%.RData: $(MODEL_BIN)
	mkdir -p results_WBScod_Age
	SCRIPT_INPUT="{ MODEL_CONFIG='$*' ; INCLUDE='commercial' ;  OUTFILE='$@' }" $(RUN_MODEL)

## Target 2
results_WBScod_Age/cod-survey-%.RData: $(MODEL_BIN)
	mkdir -p results_WBScod_Age
	SCRIPT_INPUT="{ MODEL_CONFIG='$*' ; INCLUDE='survey' ;  OUTFILE='$@' }" $(RUN_MODEL)


## Target 3
results_WBScod_Age/cod-both-%.RData: $(MODEL_BIN)
	mkdir -p results_WBScod_Age
	SCRIPT_INPUT="{ MODEL_CONFIG='$*' ; INCLUDE='both' ;  OUTFILE='$@' }" $(RUN_MODEL)

## Cleanup
clean:
	rm -f LGNB.o LGNB.so LGNB.dll
