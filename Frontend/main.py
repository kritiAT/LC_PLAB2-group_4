import pandas as pd
from group4_template.sequence_assembly import *
from group4_template.gene_finder import *
from group4_template.Blast import *
from werkzeug.utils import secure_filename
from req import UPLOAD_FOLDER
from flask import Flask, flash, render_template, request, redirect, url_for

ALLOWED_EXTENSIONS = {"fasta", "txt", "fastq"}

app = Flask(__name__)
app.secret_key = b'_5#034587Q564482c]/'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_PATH'] = 16 * 1024 * 1024


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
    import os
    file = request.files['file']
    if request.method == 'POST':
        if file.filename == '':
            flash('No file selected')
            return render_template('upload.html')
        elif file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            txtfile = open(str(os.path.join(UPLOAD_FOLDER, 'filename.txt')), "w")
            txtfile.write(filename)
            txtfile.close()

            flash('File uploaded successfully')
            result1 = get1()  # Assembled DNA
            print("Result1 done")
            obj=Translate(dna=result1, threshold=70)
            result2=''.join(obj.f_mrna)  # mRNA seq
            print("Result2 done")
            result3=''.join(obj.r_mrna)  # Reversed mRNA seq
            #obj=Translate(dna=result1)
            print("Result3 done")
            result22= obj.proteins
            print("Result4 done", result22)

            txtfile = open(str(os.path.join(UPLOAD_FOLDER, 'ORFList.txt')), "w")
            txtfile.write(str(result22))
            txtfile.close()

            final_dict = Blast_orfs(obj.proteins)

            #final_dict={'TPLALKLNLSASPLTAAKQIRIQVARLK': ['7ABG_A4','7ABF_A4', '4Y98_A'],'SFSTPASISPTVAVCSRIRPSLLPPSTRVSQVLSTTRPLIFN': ['7ABG_A4','7ABF_A4', '4Y98_A']}
            obj.proteins_table['predicted_proteins'] = obj.proteins_table['amino_acid_sequence'].apply(lambda x: final_dict[x])


            print(final_dict)
            result4=obj.proteins_table
            #result4 = pd.DataFrame.from_dict(final_dict)
            return render_template('upload.html', rs1=result1, rs2=result2,rs3=result3, tables=[result4.to_html(classes='data', header="true")])
        else:
            flash('Not allowed')
            return render_template('upload.html')
    return redirect(url_for('home'))


def get1():
    txtfile = open(str(os.path.join(UPLOAD_FOLDER, 'filename.txt')), "r")
    txtfiledata = txtfile.read()
    txtfile.close()
    fastafilepath = os.path.join(UPLOAD_FOLDER, txtfiledata)
    #obj = Assembly(sequences=str(fastafilepath))
    obj = Assembly(sequences=str(fastafilepath))
    seq = obj.assembled_sequence
    return seq



if __name__ == '__main__':
    app.run(debug=True)
