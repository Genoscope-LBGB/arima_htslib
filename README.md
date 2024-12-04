# Installation

You will need the [Meson](https://mesonbuild.com/index.html) build system to compile the code. After it has been installed
run the following commands:
```
meson setup build --prefix=$(pwd)/install -Dwrap_mode=forcefallback
meson compile -C build
meson install -C build
```

If everything went as expected, the binaries should be in `install/bin`.
