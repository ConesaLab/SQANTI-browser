FROM --platform=linux/amd64 continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Create a non-root user
RUN groupadd -r sqanti && useradd -r -g sqanti -s /bin/bash sqanti

# Copy environment file
COPY environment.yml .

# Install dependencies
RUN conda env create -f environment.yml && conda clean -afy

# Copy the rest of the application code
COPY . .

# Activate the environment and install the package
SHELL ["conda", "run", "-n", "sqanti-browser", "/bin/bash", "-c"]
RUN pip install .

# Change ownership to non-root user
USER sqanti

# Set the PATH to include the conda environment bin directory
ENV PATH /opt/conda/envs/sqanti-browser/bin:$PATH

# Set entrypoint to the installed command
ENTRYPOINT ["sqanti_browser"]

# Default command
CMD ["--help"]
