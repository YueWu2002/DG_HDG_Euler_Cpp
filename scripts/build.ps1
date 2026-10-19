param(
    [string]$Cxx = 'g++',
    [string]$Eigen = '.deps/eigen-3.4.0'
)
$ErrorActionPreference = 'Stop'
Push-Location (Join-Path $PSScriptRoot '..')
try {
    if (-not (Test-Path (Join-Path $Eigen 'Eigen/Core'))) {
        throw 'Eigen not found. Supply -Eigen <directory containing Eigen/Core>.'
    }
    New-Item -ItemType Directory -Force build/objects | Out-Null
    $flags = @('-std=c++17', '-O2', '-Wall', '-Wextra', '-Wpedantic',
               '-UNDEBUG', '-DEIGEN_RUNTIME_NO_MALLOC', '-DEIGEN_DONT_PARALLELIZE',
               '-Idg_hdg_1d/include', '-Ivendor/fastgl', '-isystem', $Eigen)
    $sources = @('src/config.cpp', 'src/physics.cpp', 'src/simple_wave.cpp',
                 'src/discretization.cpp', 'src/volume_jacobian.cpp', 'src/dg_operator.cpp',
                 'src/hybrid.cpp', 'src/time_integrator.cpp',
                 'src/linalg/block_graph.cpp', 'src/linalg/block_ilu.cpp', 'src/linalg/gmres.cpp')
    $objects = @()
    & $Cxx @flags -c vendor/fastgl/fastgl.cpp -o build/objects/fastgl.o
    if ($LASTEXITCODE -ne 0) { throw 'FastGL compilation failed' }
    $objects += 'build/objects/fastgl.o'
    foreach ($source in $sources) {
        $object = 'build/objects/' + [IO.Path]::GetFileNameWithoutExtension($source) + '.o'
        $source = 'dg_hdg_1d/' + $source
        Write-Host "Compiling $source"
        & $Cxx @flags -c $source -o $object
        if ($LASTEXITCODE -ne 0) { throw "Compilation failed: $source" }
        $objects += $object
    }
    & $Cxx @flags dg_hdg_1d/src/main.cpp @objects -o build/dg_hdg.exe
    if ($LASTEXITCODE -ne 0) { throw 'Solver link failed' }
    & $Cxx @flags -c dg_hdg_1d/check/dg_reference.cpp -o build/objects/dg_reference.o
    if ($LASTEXITCODE -ne 0) { throw 'Reference kernel compilation failed' }
    $checkObjects = $objects + @('build/objects/dg_reference.o')
    & $Cxx @flags dg_hdg_1d/check/main.cpp dg_hdg_1d/check/check_1d.cpp dg_hdg_1d/check/linearization.cpp @checkObjects -o build/check_1d.exe
    if ($LASTEXITCODE -ne 0) { throw 'Check link failed' }
    & $Cxx @flags dg_hdg_1d/check/benchmark_dg.cpp @checkObjects -o build/benchmark_dg.exe
    if ($LASTEXITCODE -ne 0) { throw 'Benchmark link failed' }
} finally {
    Pop-Location
}
