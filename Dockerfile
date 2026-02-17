FROM python:3.13-slim

ENV DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.authors="nipper@wisc.edu"

ENV DEBIAN_FRONTEND=noninteractive

COPY ./requirements.txt /app/
COPY ./Dockerfile /app/
COPY ./app.py /app/
COPY ./shared.py /app/
COPY ./calculations/ /app/calculations/
COPY ./modules/ /app/modules/

WORKDIR /app

RUN apt-get update && apt-get upgrade -y

RUN pip install --upgrade pip && pip install -r requirements.txt

EXPOSE 8000/tcp

CMD ["shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "8000"]