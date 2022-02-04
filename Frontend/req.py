import os

from pathlib import Path
from flask import Flask, flash, render_template, request, redirect, url_for

home_dir = Path.home()
PROJECT_DIR = home_dir.joinpath(".project3")
DATA_DIR = PROJECT_DIR.joinpath("data")
UPLOAD_FOLDER = os.path.join(DATA_DIR, 'uploads')
faapath = os.path.join(UPLOAD_FOLDER, 'file.faa')