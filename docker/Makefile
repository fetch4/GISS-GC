# Makefile for building, pulling, running, and pushing Docker images for GISS-GC
# development environment.
#
# The GNU compilers will be used by default. To build with Intel compilers, call
#
#   make build COMPILER=intel

# Default target
all: build

# Variables
COMPILER	?= gcc
DOCKERFILE	= Dockerfile.$(COMPILER)
NAMESPACE	= ghcr.io/fetch4
IMAGE_NAME	= giss-gc-dev-env-$(COMPILER)
TAG			= latest

# Full image name
IMAGE		= $(NAMESPACE)/$(IMAGE_NAME):$(TAG)

# Pull the Docker image from the registry
pull:
	docker pull $(IMAGE)

# Build the Docker image using the specified Dockerfile
build:
	docker build -f $(DOCKERFILE) -t $(IMAGE) .

# Run the Docker container interactively with the home directory mounted
# NOTE: This exposes your home directory to the Docker image and should be used
#       with caution.
run:
	docker run --rm -it -v ${HOME}:${HOME} $(IMAGE)

# Push the Docker image to the registry
push:
	docker push $(IMAGE)
