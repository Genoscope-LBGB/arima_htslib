# Only for Genoscope internal usage, please use the build instructions in the
# Readme file for compilation commands 

all: build-bin

build-bin:
	. /env/products/fgtools/2.0/bin/fg_bashfix && \
        module load meson gcc/12 ninja-build && \
        meson setup build --prefix=$(pwd)/install -Dwrap_mode=forcefallback && \
        meson compile -C build && \
        meson install -C build

install:
	cp -r install/* $(prefix)/