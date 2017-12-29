#!/bin/bash
docker run -ti --mount type=bind,src=/home/petiaccja/UppsalaUniversity/UU_AppliedFEM/Part3,dst=/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:latest
