default:
	@echo "Only used in make release"

release:
	make bumpversion
	git push --tags
	mkdir -p dist-old
	-mv dist/* dist-old
	python setup.py sdist
	twine upload dist/*
