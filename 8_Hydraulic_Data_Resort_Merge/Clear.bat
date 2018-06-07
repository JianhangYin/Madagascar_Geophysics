@echo off 
for /f "tokens=*" %%a in ('dir /b /ad /s D:\Paper\Data\CottonValley\resorted\^|sort /r') do rd "%%a" /q 2>nul