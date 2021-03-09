FROM python:3

WORKDIR /Simulinac
ENV DISPLAY=192.168.1.52:0
VOLUME /Simulinac

#COPY requirements.txt ./
#RUN python -m pip install --no-cache-dir -r requirements.txt
RUN python -m pip install --no-cache-dir matplotlib numpy PyQt5 PyYAML scipy
#copy . ./

CMD ["python", "simu.py"]
