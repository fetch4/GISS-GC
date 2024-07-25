all: build

DOCKERFILE 	= Dockerfile.devenv
NAMESPACE 	= ghcr.io/fetch4
IMAGE_NAME	= giss-gc-dev-env
TAG 				= latest

IMAGE				= $(NAMESPACE)/$(IMAGE_NAME):$(TAG)

pull:
	docker pull $(IMAGE)

build:
	docker build -f $(DOCKERFILE) -t $(IMAGE) --no-cache .

run:
	docker run --rm -it -v ${HOME}:${HOME} $(IMAGE)

push:
	docker push $(IMAGE)
