#!/bin/sh -eux

clang-format --version

git checkout -B formatter-work origin/fix/postgres-bingo-regression-tests-lint

clang-format -i \
  bingo/postgres/src/pg_am/pg_bingo_build.cpp \
  bingo/postgres/src/pg_am/pg_bingo_update.cpp

if ! git diff --quiet -- \
  bingo/postgres/src/pg_am/pg_bingo_build.cpp \
  bingo/postgres/src/pg_am/pg_bingo_update.cpp
do
  git config user.name "Bingo clang-format diagnostic"
  git config user.email "noreply@github.com"
  git add \
    bingo/postgres/src/pg_am/pg_bingo_build.cpp \
    bingo/postgres/src/pg_am/pg_bingo_update.cpp
  git commit -m "Apply clang-format 11.0.1 to Bingo PostgreSQL changes"
  git push origin HEAD:fix/postgres-bingo-regression-tests-lint
fi

exit 1
