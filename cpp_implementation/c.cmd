

@echo off

rem IN FLAGS:
set "main=main.cpp"
set "GPP=C:\mingw64\bin\g++.exe"

set "eigen=C:\_trabalho\codes\new_SEM\eigen-5.0.0"

rem OUT FLAGS:
set "out=sem.exe"

if exist "%out%" (
    del "%out%"
    echo file "%out%" deleted
)

echo Starting compilation
"%GPP%" "%main%" -I"%eigen%" -o "%out%"

echo Starting program
echo.

if exist "%out%" (
    "%out%"
)



