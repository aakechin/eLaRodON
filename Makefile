# Makefile for elarodon package management
include .env

# Configuration
PACKAGE_NAME := elarodon
PYTHON := python3
TWINE := twine

.PHONY: help build check deploy clean

help:
	@echo "Usage:"
	@echo "  make build    Build package"
	@echo "  make deploy   Upload to PyPI (uses ${ENV_FILE})"
	@echo "  make clean    Clean build artifacts"
	@echo ""
	@echo "Example:"
	@echo "  make deploy REPOSITORY=testpypi"

build:
	@echo "Building package..."
	${PYTHON} -m pip install --upgrade build
	${PYTHON} -m build
	@echo "\n✅ Build complete. Files:"
	@ls -lh dist/


check:
	@echo "Checking package metadata..."
	${PYTHON} -m pip install --upgrade twine
	${TWINE} check dist/*
	@echo "\n✅ Check passed. Package ready for deployment."

deploy:
	@echo "Preparing to upload to $(REPOSITORY)..."
	@test -n "$(PYPI_TOKEN)" || { echo "❌ ERROR: PYPI_TOKEN is not set"; exit 1; }
	@echo "Uploading using token: $(shell echo $(PYPI_TOKEN) | cut -c1-4)...$(shell echo $(PYPI_TOKEN) | rev | cut -c1-4 | rev)"

	@TWINE_USERNAME=__token__ TWINE_PASSWORD="$(PYPI_TOKEN)" \
		$(TWINE) upload dist/*
	
# 	@TWINE_USERNAME=__token__ TWINE_PASSWORD="$(PYPI_TOKEN)" \
# 		$(TWINE) upload dist/elarodon-0.1.5-py3-none-any.whl
	
# 	@TWINE_USERNAME=__token__ TWINE_PASSWORD="$(PYPI_TOKEN)" \
# 		$(TWINE) upload dist/elarodon-0.1.5.tar.gz
	
	
	@echo "\n✅ Upload complete! Package now available on:"
	@test "$(REPOSITORY)" = "pypi" \
		&& echo "https://pypi.org/project/$(PACKAGE_NAME)/" \
		|| echo "https://test.pypi.org/project/$(PACKAGE_NAME)/"

clean:
	@echo "Cleaning build artifacts..."
	rm -rf build/ dist/ src/elarodon/*.egg-info/
	@echo "🧹 Clean complete."