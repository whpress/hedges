.PHONY: clean

PWD= $(shell pwd)
DOCKER_IMG=bionic-hedges-python2.7
test_docker:
	docker run -v ${PWD}:/work ${DOCKER_IMG} bash -c "mkdir -p /tmp/build && cd /tmp/build && cmake /work/ && make && cd /work/ && PYTHONPATH=/work/build/src python print_module_help_files.py"

clean:
	docker run -v ${PWD}:/work ${DOCKER_IMG} bash -c "rm -r /work/build/"

build:
	docker run -v ${PWD}:/work ${DOCKER_IMG} bash -c "mkdir -p /work/build && cd /work/build && cmake /work/ && make && ctest -V"

hedges_testprogramm:
	docker run -v /home/dgerin/work/dpu_dnastorage/hedges:/work bionic-hedges-python2.7 bash -c "cd /work/ && PYTHONPATH=/work/build/src python test_program.py"

hedges_print_module_help_files:
	docker run -v /home/dgerin/work/dpu_dnastorage/hedges:/work bionic-hedges-python2.7 bash -c "cd /work/ && PYTHONPATH=/work/build/src python print_module_help_files.py"

.DEFAULT_GOAL := build
