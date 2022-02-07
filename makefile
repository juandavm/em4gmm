# Expectation Maximization for Gaussian Mixture Models.
# Copyright (C) 2012-2022 Juan Valor
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

CFLAGS := -Ofast -mtune=native -fno-math-errno -fno-trapping-math -ftree-vectorize -pipe
DFLAGS := -g -pipe
CC := gcc
ifeq (, $(shell which gcc))
    CC := clang
endif
BINDIR := /usr/bin

SRC := src
BIN := bin
DAT := dat
ALL := $(SRC)/workers.c $(SRC)/data.c $(SRC)/gmm_base.c $(SRC)/gmm_class.c $(SRC)/gmm_train.c $(SRC)/gmm_merge.c

all:
	$(CC) $(CFLAGS) $(ALL) $(SRC)/class.c -o $(BIN)/gmmclass -pthread -lz -lm
	$(CC) $(CFLAGS) $(ALL) $(SRC)/train.c -o $(BIN)/gmmtrain -pthread -lz -lm
	chmod +x $(BIN)/gmmclass
	chmod +x $(BIN)/gmmtrain
	strip $(BIN)/gmmtrain
	strip $(BIN)/gmmclass

debug:
	$(CC) $(DFLAGS) $(ALL) $(SRC)/class.c -o $(BIN)/gmmclass -pthread -lz -lm
	$(CC) $(DFLAGS) $(ALL) $(SRC)/train.c -o $(BIN)/gmmtrain -pthread -lz -lm
	chmod +x $(BIN)/gmmclass
	chmod +x $(BIN)/gmmtrain

install:
	cp -f $(BIN)/gmmclass $(BINDIR)/gmmclass
	cp -f $(BIN)/gmmtrain $(BINDIR)/gmmtrain
	cp -f $(BIN)/datagen.py $(BINDIR)/datagen.py
	chmod +x $(BINDIR)/gmmclass
	chmod +x $(BINDIR)/gmmtrain
	chmod +x $(BINDIR)/datagen.py

clean:
	rm -f $(BIN)/gmmtrain
	rm -f $(BIN)/gmmtrain
