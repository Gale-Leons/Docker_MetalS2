FROM python:3.12-slim 

RUN pip install --upgrade pip

# Installa fork patchata di Biopython + dipendenze
RUN pip install git+https://github.com/Gale-Leons/biopython@master
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copia gli script del progetto
COPY scripts/ app/scripts/

WORKDIR /app
