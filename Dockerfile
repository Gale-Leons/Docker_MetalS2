FROM python:3.12-slim 

# Install dependecies for building Python packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libffi-dev \
    python3-dev \
    git \
    procps \
    && rm -rf /var/lib/apt/lists/*


ARG USER_ID=1000
ARG GROUP_ID=1000
RUN groupadd -g ${GROUP_ID} appuser && \
    useradd -m -u ${USER_ID} -g appuser appuser


RUN mkdir -p /app && chown -R appuser:appuser /app
RUN mkdir -p /tmp && chmod 777 /tmp

USER appuser

WORKDIR /app
COPY --chown=appuser:appuser . .

# Install your fork of Biopython
RUN pip3 install --upgrade pip setuptools wheel
RUN pip3 install git+https://github.com/Gale-Leons/biopython@master

# Install other dependecies
COPY requirements.txt .
RUN pip3 install -r requirements.txt

# Default help command
CMD ["python3", "scripts/metals2.py", "--h"]
