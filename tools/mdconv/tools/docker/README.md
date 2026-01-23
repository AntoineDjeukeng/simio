# Docker Dev Environment

This container provides a clean Ubuntu environment with GROMACS and LAMMPS installed.

## Build the image

```bash
docker compose -f tools/docker/compose.yml build
```

## Enter the container

```bash
docker compose -f tools/docker/compose.yml run --rm mdconv-dev bash
```

## Build and test inside Docker

```bash
make
make test
```
