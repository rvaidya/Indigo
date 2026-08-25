#!/bin/sh -eux

clang-format --version
clang-format -i \
  bingo/postgres/src/pg_am/pg_bingo_build.cpp \
  bingo/postgres/src/pg_am/pg_bingo_update.cpp

git diff -- \
  bingo/postgres/src/pg_am/pg_bingo_build.cpp \
  bingo/postgres/src/pg_am/pg_bingo_update.cpp

exit 1
