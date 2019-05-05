LIB_DIR = lib

default: TDDispl

TDDispl: setup.py TDDispl.pyx $(LIB_DIR)/libexamples.a
	python3 setup.py build_ext --inplace && rm -f TDDispl.c && rm -Rf build

$(LIB_DIR)/libexamples.a:
	make -C $(LIB_DIR) libexamples.a

clean:
	rm *.so
