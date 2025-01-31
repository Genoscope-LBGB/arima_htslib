# Only for Genoscope internal usage, please use the build instructions in the
# Readme file for compilation commands 

all: build-bin

build-bin:
	. /env/products/fgtools/2.0/bin/fg_bashfix && \
        module load meson ninja-build cmake/3.19.1 && \
        meson setup build --prefix=$(PWD)/install.d -Dwrap_mode=forcefallback && \
        meson compile -C build && \
        meson install -C build

install:
	cp -r install.d/* $(PREFIX)/