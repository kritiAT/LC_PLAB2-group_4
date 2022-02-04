import os

from task1 import *
from task2 import *
from readingfasta import fasta_list
from req import UPLOAD_FOLDER, faapath
from flask import Flask, flash, render_template, request, redirect, url_for


ALLOWED_EXTENSIONS = {"fasta", "fna", "ffn", "faa", "frn", "fa"}

app = Flask(__name__)
app.secret_key = b'_5#034587Q564482c]/'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_PATH'] = 16 * 1024 * 1024
dna_sequence = fasta_list(faapath)
dna_sequence = dna_sequence[0][1]


@app.route("/")
def home():
    return render_template('index.html')


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/upload')
def upload():
    return render_template('upload.html')


@app.route('/uploader', methods=['GET', 'POST'])



def upload_file():

    file = request.files['file']

    if request.method == 'POST':

        if file.filename == '':
            flash('No file selected')
            return render_template('upload.html')
        elif file and allowed_file(file.filename):
            filename = "file.faa"
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            flash('File uploaded successfully')
            result1 = get1()
            result2 = list(AA_prediction(result1, 0))


            return render_template('upload.html', rs1=result1, rs2=result2)
        else:
            flash('Not allowed')
            return render_template('upload.html')

        return redirect(url_for('home'))

def get1():
    exec(open("task1.py").read())
    main_seq=seq_assembly(sequences=dna_kmers_random, score_mat=scoring_matrix)
    return main_seq




if __name__ == '__main__':
    app.run(debug=True)
