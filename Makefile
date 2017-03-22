all: dir
	@cd build && cmake ..
dir:
	@rm -rf build
	@mkdir build

clean:
	@git clean -xdf
