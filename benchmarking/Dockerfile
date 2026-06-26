FROM mambaorg/micromamba:2-ubuntu22.04

# x3cflux 3.0.1 requires Linux x86_64 (manylinux_2_27+). This image satisfies that.
# On Mac/Windows use Docker Desktop — x3cflux cannot be installed natively there.

LABEL org.opencontainers.image.title="mfa-benchmarking" \
      org.opencontainers.image.description="13C MFA: freeflux 0.3.7, influx-si 7.4.3, x3cflux 3.0.1"

USER root

# Step 1: conda-forge packages for everything that has a native build (cyipopt bundles IPOPT)
RUN micromamba install -y -n base -c conda-forge \
        python=3.10 \
        cyipopt \
        numpy">=1.20,<1.23" \
        scipy">=1.5,<1.9" \
        "pandas>=1.0,<2.0" \
        "matplotlib>=3.2,<=3.7" \
        "pyomo>=6.2,<=6.5" \
        sympy \
        openpyxl \
        xlrd \
        pillow \
        jupyterlab \
        ipykernel \
    && micromamba clean --all --yes

# Step 2: pip for PyPI-only packages (x3cflux, freeflux, influx-si, isocor, hopsy, etc.)
RUN /opt/conda/bin/pip install --no-cache-dir \
        freeflux==0.3.7 \
        influx-si==7.4.3 \
        x3cflux==3.0.1 \
        isocor==2.2.3 \
        hopsy \
        cvxpy \
        arviz

WORKDIR /workspace

EXPOSE 8888

# Mount notebooks/data at /workspace:
#   docker run -p 8888:8888 -v $(pwd)/data:/workspace/data mfa-benchmarking
CMD ["/opt/conda/bin/jupyter", "lab", \
     "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", \
     "--ServerApp.token=''", "--ServerApp.password=''"]
