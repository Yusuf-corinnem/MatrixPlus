FROM ubuntu
COPY . /project
RUN apt update && apt install gcc make vim check valgrind g++ libgtest-dev cmake -y
RUN cd /usr/src/googletest && cmake CMakeLists.txt && make