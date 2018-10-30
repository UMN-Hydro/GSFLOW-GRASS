FROM ubuntu:18.10
ENV DEBIAN_FRONTEND noninteractive
USER root

RUN useradd -m -s /bin/bash mezcla

RUN apt update && apt install tzdata

ADD install_gsflow_grass /install_gsflow_grass

RUN apt update && \
    chmod +x /install_gsflow_grass && /install_gsflow_grass

ADD install_grass_package /install_grass_package
RUN chmod +x /install_grass_package && /install_grass_package

#ADD install_binaries_grass /install_binaries_grass
#RUN chmod +x /install_binaries_grass && /install_binaries_grass

ADD install_kit /home/mezcla/install_kit
RUN chmod +x /home/mezcla/install_kit

RUN apt install -y libgl1-mesa-glx libgl1-mesa-dri
USER mezcla
WORKDIR /home/mezcla
RUN ./install_kit


ENV TERM xterm-256color
CMD sh -c "echo You need to run grass from a shell created by a docker exec invocation; sleep 5d"
