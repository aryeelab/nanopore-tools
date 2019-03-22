# Oxford nanopore tools

## Installing dependencies

- Install Homebrew

https://brew.sh

- Install cromwell

```
brew install cromwell
```

- Install Docker

https://hub.docker.com/editions/community/docker-ce-desktop-mac

## Getting started

### Quick Start

```
cromwell run -i test_data/test-run-1.json preprocess_flowcell.wdl 
```

You can include the path to a custom Cromwell config file by something like:

```
JAVA_OPTS="-Dconfig.file=/Users/maryee/cromwell/cromwell.conf"
cromwell run -i test_data/test-run-1.json preprocess_flowcell.wdl 
```
(See https://cromwell.readthedocs.io/en/latest/tutorials/ConfigurationFiles/)

This can be used to specify a cache database, for example:

```
docker run -p 3306:3306 --name cromwell_mysql -e MYSQL_ROOT_PASSWORD=pass123 -e MYSQL_DATABASE=cromwell -e MYSQL_USER=cromwell -e MYSQL_PASSWORD=pass123 -d mysql/mysql-server:5.5

```
(See https://cromwell.readthedocs.io/en/latest/tutorials/PersistentServer/)

### Testing within a docker image

e.g.
```
docker run --rm -it aryeelab/guppy
```

## Visualizing the workflow graph

You can use `womtool` (part of cromwell) to output a workflow graph in `.dot` format:

```
womtool graph preprocess_flowcell.wdl > preprocess_flowcell.dot
```

This graph can be edited if necessary (such as to label edges with inputs/outputs), and then visualized with graphviz (`brew install graphviz`):

```
dot preprocess_flowcell.dot -Tpng -o preprocess_flowcell.png
```

The workflow graph below is produced in this way.


## Preprocess Flowcell workflow

![alt text](preprocess_flowcell.png "preprocess_flowcell.wdl DAG")

