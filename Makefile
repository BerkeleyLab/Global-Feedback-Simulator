all:
	sh doswig.sh

clean:
	rm -f *.o *_wrap.c *.so *.pyc *.so linac.py filter.py
