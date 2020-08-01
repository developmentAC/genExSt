FROM python:3
# see notes from ref https://discuss.streamlit.io/t/how-to-use-streamlit-in-docker/1067/7
# RUN apt-get update

RUN \
    apt-get update &&\
    apt-get install -y git &&\
    apt-get install -y htop &&\
    apt-get install -y vim &&\
		apt-get install -y python3 &&\
		apt-get install -y python3-pip &&\
		pip install streamlit &&\
		pip install pyvis &&\
#		pip install spacy &&\
		pip install pytest &&\
		pip install plotly_express &&\
		pip install pip install jsonpickle &&\
		pip install scikit-learn &&\
		pip install chart_studio &&\
		pip install matplotlib

		RUN useradd mydevi
		RUN mkdir /home/mydevi
		RUN export HOME=/home/mydevi

WORKDIR /home/mydevi
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN mkdir -p /root/.streamlit

RUN bash -c 'echo -e "\
[general]\n\
email = \"\"\n\
" > /root/.streamlit/credentials.toml'

RUN bash -c 'echo -e "\
[server]\n\
enableCORS = false\n\
" > /root/.streamlit/config.toml'

EXPOSE 8501

CMD bash


# Summary of commands
# Build: sudo docker build -t mydevi .
# Mount local directory and run container:
# sudo docker run -it -p 8501:8501 --mount type=bind,source=$PWD,target=/home/mydevi mydevi
