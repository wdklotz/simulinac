FROM python:3

WORKDIR /Simulinac
ENV DISPLAY=192.168.1.52:0
VOLUME /Simulinac

RUN python -m pip install --no-cache-dir matplotlib PyYAML scipy 

CMD ["python", "simu.py"]
