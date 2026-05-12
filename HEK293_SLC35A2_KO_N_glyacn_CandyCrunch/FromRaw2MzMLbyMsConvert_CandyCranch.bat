@echo off

for /D %%f in (*.raw) do (
    echo Processing: %%f
    msconvert "%%f" --filter "peakPicking true 1-" --32 --mzXML -o . --outfile "%%~nf.mzXML" -v
)

echo All done!
pause