@echo off
mkdir test 2>nul >nul
pushd test
if not exist "rt.h" (
	curl -LJO https://raw.githubusercontent.com/leok7v/ui/main/single_file_lib/rt/rt.h
)
if not exist "ui.h" (
	curl -LJO https://raw.githubusercontent.com/leok7v/ui/main/single_file_lib/ui/ui.h
)
if not exist "sqlite3.c" (
	curl -LJO https://raw.githubusercontent.com/jmscreation/libsqlite3/main/src/sqlite3.c
)
popd

