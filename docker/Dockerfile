# DFM CAPE engine: cape_request.json → cape_result.json (no COM).
#
# Build from this folder (the Aspen machine payload):
#   docker build -f docker/Dockerfile -t dfm-methanation-cape:v1 .
#
# Run:
#   docker run --rm -v "$PWD:/data" dfm-methanation-cape:v1 \
#     /data/cape_request.json -o /data/cape_result.json
#
# Host dolfinx is 0.10.0.post3; pin the matching end-user image, not :stable (0.11).

ARG DOLFINX_IMAGE=dolfinx/dolfinx:v0.10.0-r1

FROM ${DOLFINX_IMAGE} AS cape-engine

LABEL org.opencontainers.image.title="DFM methanation CAPE engine" \
      org.opencontainers.image.description="CAPE-OPEN 1.1 stand-in: cape_request_v1 in, cape_result_v1 out"

WORKDIR /opt/dfm

ENV PYTHONPATH=/opt/dfm:$PYTHONPATH \
    PYTHONUNBUFFERED=1 \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

# Whitelist only. Not copied: fermentation, system_c, plots/web/Optuna,
# COMSOL extracts, lab xlsx, MPH, bulky results/.
COPY adr/__init__.py adr/__init__.py
COPY adr/cape/ adr/cape/
COPY adr/controller/__init__.py adr/controller/simulation.py adr/controller/
COPY adr/factories/__init__.py adr/factories/base.py \
     adr/factories/methanation_dfm.py adr/factories/registry.py adr/factories/
COPY adr/model/ adr/model/
COPY adr/numerical/__init__.py adr/numerical/assembler.py \
     adr/numerical/dfm_runtime.py adr/numerical/solver.py adr/numerical/
COPY adr/view/__init__.py adr/view/export.py adr/view/
COPY fenicsx/__init__.py fenicsx/balance_diagnostics.py \
     fenicsx/boundary_conditions.py fenicsx/cycle_scheduler.py \
     fenicsx/darcy_mass_balance.py fenicsx/diluted_species.py \
     fenicsx/inlet_sopdt.py fenicsx/mesh_cylinder.py \
     fenicsx/properties.py fenicsx/run_export.py \
     fenicsx/wilke_transport.py fenicsx/
COPY dfm_config.py geometry_config.py kinetics.py ./
COPY config/dfm_config.pow_v07_10_60_opt-1.json \
     config/cape_request.example.json config/

ENTRYPOINT ["python3", "-m", "adr.cape"]
CMD ["--help"]
