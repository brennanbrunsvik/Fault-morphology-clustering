LIB_DIR = lib

default: TDStrain

TDStrain: setup.py TDStrain.pyx $(LIB_DIR)/libexamples.a
	python3 setup.py build_ext --inplace && rm -f TDStrain.c && rm -Rf build

$(LIB_DIR)/libexamples.a:
	make -C $(LIB_DIR) libexamples.a

clean:
	rm *.so
