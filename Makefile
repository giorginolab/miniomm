default:
	@echo "Only used in make release"

release:
	bumpversion patch
	git push --tags
	mkdir -p dist-old
	-mv dist/* dist-old
	python -m build
	twine upload dist/*
