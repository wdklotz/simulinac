FROM python:3

WORKDIR ${PWD}

COPY requirements.txt ./
RUN python -m pip install --no-cache-dir -r requirements.txt

copy . .

CMD ["python", "simu.py"]
