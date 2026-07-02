@echo off

:: Write the config file
(
echo filter=peakPicking vendor snr=1 peakSpace=0.1 msLevel=1-
echo filter=turbocharger halfIsoWidth=2
echo filter=titleMaker ^<RunId^>.^<ScanNumber^>.^<ScanNumber^>.^<ChargeState^> File:"^<SourcePath^>", NativeID:"^<Id^>"
echo filter=threshold tic-relative 0.001 most-intense
) > msconvert_config.txt

for /D %%f in (*.raw) do (
    echo Processing: %%f
    msconvert --zlib --64 -c msconvert_config.txt -o . --outfile "%%~nf_denoised_by_msconvert.mzML" -v "%%f"
)

del msconvert_config.txt
echo All done!
pause







::--filter "threshold count 100 least-intense


::MS2Denoise [<peaks_in_window> [<window_width_Da> [multicharge_fragment_relaxation]]]

