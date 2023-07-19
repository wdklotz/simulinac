# syntax=docker/dockerfile:1
FROM python:3.11
WORKDIR /src
COPY . .
ENV DISPLAY=192.168.1.52:0
#VOLUME /home/wdklotz/SIMULINAC_Copy

RUN python3 -m pip install --no-cache-dir matplotlib PyYAML scipy 

CMD ["python", "simu.py"]
