# Getting started

:warning: **This app is still under development!** :warning:

## Run locally

### Server

Create a local environment with conda and install server side dependencies with pip from `src/server/requirements.txt`:

```bash
conda create -n moltools python=3.10
conda activate moltools
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