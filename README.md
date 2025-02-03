[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Maintainer](https://img.shields.io/badge/Maintainer-davidmeijer-blue)](https://github.com/davidmeijer)
[![Generic badge](https://img.shields.io/badge/Version-alpha-green.svg)](https://shields.io/)

# MolTools

<img src="./logo.png" alt="logo" width="100">

MolTools is a web-based toolbelt for visualizing and analyzing natural product compounds.

Visit the MolTools toolbelt [here](https://moltools.bioinformatics.nl/).

## Run locally for development

### Server

Create a local environment and install server side dependencies with pip from `src/server/requirements.txt`:

```bash
pip install -r src/server/requirements.txt
```

### Client

First install NPM package manager and Node.js on your device.

Then install client side dependencies with NPM from `src/client/package.json`:

```bash
cd src/client
npm install
```

### Run

Run the server in one terminal:

```bash
python3 ./app/server/api.py
```

Run the client in another terminal:

```bash
cd src/client
npm start
```

Visit `http://localhost:3000/` to view the app.

## Run with Docker

Run the following script to build and runt he app in a Docker container:

```bash
docker-compose -p moltools up --build --force-recreate --remove-orphans -d
```

The app will be available at `https://localhost:4001/`.