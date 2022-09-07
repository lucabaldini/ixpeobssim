docs: html pdf

package:
	cd docs; ./package_docs.py

html:
	cd docs; python macro/datafmt.py; make html

htmlall:
	cd docs; python macro/datafmt.py; make htmlall

pdf:
	cd docs; python macro/datafmt.py; make latex; make latexpdf

test:
	cd tests; make

clean:
	./tools/cleanup.py
	rm -rf localinstall

localinstall:
	mkdir localinstall; python setup.py install --prefix localinstall

