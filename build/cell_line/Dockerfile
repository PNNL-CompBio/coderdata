FROM r-base:4.0.3
#ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update --allow-insecure-repositories
RUN apt-get install -y --allow-unauthenticated build-essential libpq-dev python3.10 python3-pip python3-setuptools python3-dev python3-venv

RUN python3 -m venv /opt/venv
RUN /opt/venv/bin/pip3 install --upgrade pip

ENV PYTHONPATH "${PYTHONPATH}:/app"
WORKDIR /app

ADD requirements.txt .
ADD requirements.r .
ADD *.R .
ADD *.py .

# installing python libraries
RUN /opt/venv/bin/pip3 install -r requirements.txt

# installing r libraries
RUN Rscript requirements.r

