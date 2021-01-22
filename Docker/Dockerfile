FROM comics/centos8

RUN dnf upgrade -y
RUN dnf -y groupinstall "Development Tools"
RUN dnf -y install openssl-devel bzip2-devel libffi-devel
RUN dnf -y install python3.8

RUN pip3.8 install numpy configparser matplotlib scipy

RUN dnf -y install python38-tkinter



