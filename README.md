# Oxford nanopore tools

## Install dependencies

- Install Homebrew

https://brew.sh

- Install cromwell

```
brew install cromwell
```

- Install Docker

https://hub.docker.com/editions/community/docker-ce-desktop-mac


## Quick Start

```
cromwell run -i test_data/test-run-1.json demultiplex_and_basecall.wdl 
```


## Testing within a docker image

e.g.
```
docker run --rm -it aryeelab/nanopore_albacore
```