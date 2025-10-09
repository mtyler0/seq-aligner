FROM python:3.14.0rc3-slim-trixie
WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir -r requirements.txt
RUN mkdir -p /src/data/uploaded_files
RUN mkdir -p /src/resources
EXPOSE 5000
CMD ["python", "src/app.py"]