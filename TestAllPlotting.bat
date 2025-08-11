@echo off

:: Batch script to run through adjoint test cases

setlocal enabledelayedexpansion

REM Define the parameters
set tests=tanh_perflux_ke-2.h5 randflux2_ke-2.h5
$middleFlux.h5 
$set mins=Saddles SimpleStep
set mins=Saddles 
$SimpleStep
set stepFraction=0.8 0.6 0.4 0.2

REM Loop through the parameters and run the Python script
for %%t in (%tests%) do (
    echo Setting minimization style

    
    for %%s in (%mins%) do (

        if "%%s"=="Saddles" (
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







