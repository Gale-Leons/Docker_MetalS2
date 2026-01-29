FROM python:3.12-slim 

# Install dependecies for building Python packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libffi-dev \
    python3-dev \
    git \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /app
COPY . .

# Install your fork of Biopython
RUN pip install --upgrade pip setuptools wheel
RUN pip install git+https://github.com/Gale-Leons/biopython@master

# Install other dependecies
COPY requirements.txt .
RUN pip install -r requirements.txt

# Default help command
CMD ["python", "scripts/metals2.py", "--h"]
