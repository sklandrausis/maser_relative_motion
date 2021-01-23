FROM comics/centos8

RUN dnf upgrade -y
RUN dnf -y groupinstall "Development Tools"
RUN dnf -y install python3.8
RUN dnf -y install qt5-qtbase qt5-qtbase-devel
RUN dnf -y install mesa-libGL-devel

RUN pip3.8 install numpy configparser matplotlib scipy pyqt5==5.14

RUN dnf -y install python38-tkinter