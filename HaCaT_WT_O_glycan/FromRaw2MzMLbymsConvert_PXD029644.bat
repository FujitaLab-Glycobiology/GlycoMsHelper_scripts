@echo off

:: Write the config file
(
echo filter=peakPicking vendor snr=1 peakSpace=0.1 msLevel=1-
echo filter=turbocharger halfIsoWidth=2
echo filter=titleMaker ^<RunId^>.^<ScanNumber^>.^<ScanNumber^>.^<ChargeState^> File:"^<SourcePath^>", NativeID:"^<Id^>"
) > msconvert_config.txt

for %%f in (*.raw) do (
    echo Processing: %%f
    msconvert --zlib --64 -c msconvert_config.txt -o . --outfile "%%~nf.mzML" -v "%%f"
)

del msconvert_config.txt
echo All done!
pause