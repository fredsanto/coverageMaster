FROM python:3.7-slim

RUN pip3 install numpy==1.16.2 sympy==1.0 PyWavelets==1.0.3 matplotlib==2.2.3 scipy==1.2.1 more_itertools==8.8.0 pandas==1.0.0

COPY . .

ENTRYPOINT ["python","coverageMaster.py"]
