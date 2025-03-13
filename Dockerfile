FROM python:3.9 AS base

# install bedtools
RUN apt-get update
RUN apt-get install bedtools

# install python packages
COPY . /app
RUN pip install --upgrade pip && \
pip install --trusted-host pypi.python.org -r /app/requirements.txt

# set environment variables
ENV PATH=$PATH:/app
RUN chmod a+x /app/WGS.py
RUN chmod a+x /app/SOS.py

WORKDIR /mnt