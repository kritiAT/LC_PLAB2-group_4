# syntax=docker/dockerfile:1
FROM python:3.10
LABEL maintainer="Sarah Imani Moghaddam <s0saiman@uni-bonn.de>, Elena Xerxa <elena.xerxa@gmail.com>, Kriti Amin <kritiamin6461@qmail.com>, Faiza Khurshidi <faiza.khurshid@uni-bonn.de>"
WORKDIR /app
COPY . /app
RUN pip3 install group4_template/
ENV FLASK_PORT 5000
EXPOSE ${FLASK_PORT}
#docker build . -t group4:latest
#docker run --name group4_run1 -p 5000:5000 -d group4:latest
ENTRYPOINT [ "python",  "main.py"  ]