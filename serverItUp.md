# Summary of commands to create a Docker container for this project.
# Build:
sudo docker build -t mydevi .
# Mount local directory and run container:
sudo docker run -it -p 8501:8501 --mount type=bind,source=$PWD,target=/home/mydevi mydevi
