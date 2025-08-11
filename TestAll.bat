@echo off

:: Batch script to run through adjoint test cases

setlocal enabledelayedexpansion

REM Define the parameters
set tests=perflux2_ke-2.h5  
$ extraDiff2_ke-2+2e-3.h5 randflux2_ke-2.h5 
set mins=Saddles
$GSS SimpleStep
set stepFraction=1.0

REM Loop through the parameters and run the Python script
for %%t in (%tests%) do (
    echo Setting minimization style

    
    for %%s in (%mins%) do (

        if "%%s" NEQ "SimpleStep" (
            echo Running script with %%s for %%t
            python Master.py %%t %%s
        ) else (
            
            for %%f in (%stepFraction%) do (
                echo Running script with %%s at %%f for %%t
                python Master.py %%t %%s %%f
            )
        )
    )
)
pause
endlocal







