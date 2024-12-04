FROM python:3.13-slim

COPY ./requirements.txt /
COPY ./Dockerfile /
COPY ./app.py /

RUN apt-get update && apt-get upgrade 

RUN pip install --upgrade pip && pip install -r requirements.txt

EXPOSE 8000

CMD ["shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "8000"]