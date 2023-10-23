FROM python:3.11
EXPOSE 8500
COPY data/ ./ \
     config.json ./
RUN apt-get update && apt-get upgrade -y \
  && python -m pip install --upgrade pip \
  && pip install .

ENTRYPOINT ["streamlit", "run", "/timering/app.py", "--server.port=8500", "--server.address=0.0.0.0", "--", "--config", "config.json"]
