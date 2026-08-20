FROM python:3.10-slim

WORKDIR /app

COPY setup.py pyproject.toml ./
COPY coveragemaster/ ./coveragemaster/

RUN pip install --no-cache-dir -e .

ENTRYPOINT ["coveragemaster"]
