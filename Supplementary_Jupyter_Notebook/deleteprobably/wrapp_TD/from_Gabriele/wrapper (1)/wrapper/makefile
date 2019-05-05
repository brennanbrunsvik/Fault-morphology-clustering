LIB_DIR = lib

default: pywrapper

pywrapper: setup.py pywrapper.pyx $(LIB_DIR)/libwrapper.a
	python setup.py build_ext --inplace && rm -f pywrapper.c && rm -Rf build

$(LIB_DIR)/libwrapper.a:
	make -C $(LIB_DIR) libwrapper.a

clean:
	rm *.so
