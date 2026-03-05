FROM python:3.12-slim
WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir -r requirements.txt
RUN mkdir -p app/src/data/uploaded_files
RUN mkdir -p app/src/resources
EXPOSE 5000
CMD ["python", "src/app.py"]