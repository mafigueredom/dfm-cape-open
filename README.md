# DFM Methanation Reactor — Aspen Plus payload

CAPE-OPEN 1.1 unit for 64-bit Aspen Plus. This folder is the whole module: the Docker engine and the COM DLL. It does not depend on the rest of the doctoral-research tree.

Palette name: **DFM Methanation Reactor**  
ProgID: `PhD.DFM.Methanation.1`  
Image: `dfm-methanation-cape:v1`

## On the Aspen PC

Install Docker Desktop (leave it running), the .NET SDK, and .NET Framework 4.8. Aspen Plus must be 64-bit.

From this folder:

```bat
build.bat
```

Then open an **elevated** command prompt in this folder and run:

```bat
com\install.bat
```

Restart Aspen Plus. Place **DFM Methanation Reactor** from the CAPE-OPEN / user-model palette.

Connect `Feed_ads`, `Feed_purge`, `Feed_rxn`, and `Product`. Connect `Feed_purge2` only if `t_purge2` is greater than 0. The property package must include CO2 `124-38-9`, H2 `1333-74-0`, CH4 `74-82-8`, H2O `7732-18-5`, and N2 `7727-37-9`.

`Calculate()` runs `docker` in the same Windows session, so Docker Desktop has to be running. The first build pulls `dolfinx/dolfinx:v0.10.0-r1`.

Unregister with `com\uninstall.bat` (elevated).
