FROM python:3.12-slim 

# Install dependecies for building Python packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libffi-dev \
    python3-dev \
    git \
    procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . . 

# Install your fork of Biopython
RUN pip3 install --upgrade pip setuptools wheel
RUN pip3 install git+https://github.com/Gale-Leons/biopython@master
COPY requirements.txt app/requirements.txt
RUN pip3 install --no-cache-dir -r /app/requirements.txt

# verifica 
CMD ["python3", "-c", "import", "numpy;", "print(numpy.__file__)"]

# Cambia owner
RUN chmod -R 777 /app

RUN mkdir -p /tmp && chmod 777 /tmp

# Default help command
CMD ["python3", "scripts/metals2.py", "--h"]
